// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"

#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Matrix.hpp"
#include "DataStructures/Transpose.hpp"
#include "DataStructures/Variables.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "NumericalAlgorithms/Spectral/Spectral.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/Spherepack.hpp"
#include "NumericalAlgorithms/SphericalHarmonics/SpherepackCache.hpp"
#include "Utilities/Algorithm.hpp"
#include "Utilities/Blas.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeArray.hpp"
#include "Utilities/MemoryHelpers.hpp"
#include "Utilities/StdArrayHelpers.hpp"

namespace partial_derivatives_detail {
template <size_t Dim, typename VariableTags, typename DerivativeTags>
struct LogicalImpl;

// This routine has been optimized to perform really well. The following
// describes what optimizations were made.
//
// - The `partial_derivatives` functions below have an overload where the
//   logical derivatives may be passed in instead of being computed. In the
//   overloads where the logical derivatives are not passed in they must be
//   computed. However, it is more efficient to allocate the memory for the
//   logical partial derivatives with respect to each coordinate at once. This
//   requires the `partial_derivatives_impl` to accept raw pointers to doubles
//   for the logical derivatives so it can be used for all overloads.
//
// - The resultant Variables `du` is a not_null pointer so that mutating compute
//   items can be supported.
//
// - The storage indices into the inverse Jacobian are precomputed to avoid
//   having to recompute them for each tensor component of `u`.
//
// - The DataVectors lhs and logical_du are non-owning DataVectors to be able to
//   plug into the optimized expression templates. This requires a `const_cast`
//   even though we will never change the `double*`.
//
// - Loop over every Tensor component in the variables by incrementing a raw
//   pointer to the contiguous data (vs. looping over each Tensor in the
//   variables with a tmpl::for_each then iterating over each component of this
//   Tensor).
//
// - We factor out the `logical_deriv_index == 0` case so that we do not need to
//   zero the memory in `du` before the computation.
template <typename ResultTags, size_t Dim, typename DerivativeFrame,
          typename ValueType = typename Variables<ResultTags>::value_type,
          typename VectorType = typename Variables<ResultTags>::vector_type>
void partial_derivatives_impl(
    const gsl::not_null<Variables<ResultTags>*> du,
    const std::array<const ValueType*, Dim>& logical_partial_derivatives_of_u,
    const size_t number_of_independent_components,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>& inverse_jacobian) {
  ValueType* pdu = du->data();
  const size_t num_grid_points = du->number_of_grid_points();
  VectorType lhs{};
  VectorType logical_du{};

  std::array<std::array<size_t, Dim>, Dim> indices{};
  for (size_t deriv_index = 0; deriv_index < Dim; ++deriv_index) {
    for (size_t d = 0; d < Dim; ++d) {
      gsl::at(gsl::at(indices, d), deriv_index) =
          InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>::get_storage_index(d, deriv_index);
    }
  }

  for (size_t component_index = 0;
       component_index < number_of_independent_components; ++component_index) {
    for (size_t deriv_index = 0; deriv_index < Dim; ++deriv_index) {
      lhs.set_data_ref(pdu, num_grid_points);
      // clang-tidy: const cast is fine since we won't modify the data and we
      // need it to easily hook into the expression templates.
      logical_du.set_data_ref(
          const_cast<ValueType*>(                              // NOLINT
              gsl::at(logical_partial_derivatives_of_u, 0)) +  // NOLINT
              component_index * num_grid_points,
          num_grid_points);
      lhs = (*(inverse_jacobian.begin() + gsl::at(indices[0], deriv_index))) *
            logical_du;
      for (size_t logical_deriv_index = 1; logical_deriv_index < Dim;
           ++logical_deriv_index) {
        // clang-tidy: const cast is fine since we won't modify the data and we
        // need it to easily hook into the expression templates.
        logical_du.set_data_ref(const_cast<ValueType*>(  // NOLINT
                                    gsl::at(logical_partial_derivatives_of_u,
                                            logical_deriv_index)) +  // NOLINT
                                    component_index * num_grid_points,
                                num_grid_points);
        lhs +=
            (*(inverse_jacobian.begin() +
               gsl::at(gsl::at(indices, logical_deriv_index), deriv_index))) *
            logical_du;
      }
      // clang-tidy: no pointer arithmetic
      pdu += num_grid_points;  // NOLINT
    }
  }
}
}  // namespace partial_derivatives_detail

template <typename DerivativeTags, typename VariableTags, size_t Dim>
void logical_partial_derivatives(
    const gsl::not_null<std::array<Variables<DerivativeTags>, Dim>*>
        logical_partial_derivatives_of_u,
    const Variables<VariableTags>& u, const Mesh<Dim>& mesh) {
  using ValueType = typename Variables<VariableTags>::value_type;
  if (UNLIKELY((*logical_partial_derivatives_of_u)[0].number_of_grid_points() !=
               u.number_of_grid_points())) {
    for (auto& deriv : *logical_partial_derivatives_of_u) {
      deriv.initialize(u.number_of_grid_points());
    }
  }
  std::array<ValueType*, Dim> deriv_pointers{};
  for (size_t i = 0; i < Dim; ++i) {
    gsl::at(deriv_pointers, i) =
        gsl::at(*logical_partial_derivatives_of_u, i).data();
  }
  if constexpr (Dim == 1) {
    Variables<DerivativeTags>* temp = nullptr;
    partial_derivatives_detail::LogicalImpl<Dim, VariableTags, DerivativeTags>::
        apply(make_not_null(&deriv_pointers), temp, temp, u, mesh);
    return;
  } else {
    auto buffer = cpp20::make_unique_for_overwrite<ValueType[]>(
        2 * u.number_of_grid_points() *
        Variables<DerivativeTags>::number_of_independent_components);
    Variables<DerivativeTags> temp0(
        &buffer[0],
        u.number_of_grid_points() *
            Variables<DerivativeTags>::number_of_independent_components);
    Variables<DerivativeTags> temp1(
        &buffer[u.number_of_grid_points() *
                Variables<DerivativeTags>::number_of_independent_components],
        u.number_of_grid_points() *
            Variables<DerivativeTags>::number_of_independent_components);
    partial_derivatives_detail::LogicalImpl<Dim, VariableTags, DerivativeTags>::
        apply(make_not_null(&deriv_pointers), &temp0, &temp1, u, mesh);
  }
}

template <typename DerivativeTags, typename VariableTags, size_t Dim>
std::array<Variables<DerivativeTags>, Dim> logical_partial_derivatives(
    const Variables<VariableTags>& u, const Mesh<Dim>& mesh) {
  auto logical_partial_derivatives_of_u =
      make_array<Dim>(Variables<DerivativeTags>(u.number_of_grid_points()));
  logical_partial_derivatives<DerivativeTags>(
      make_not_null(&logical_partial_derivatives_of_u), u, mesh);
  return logical_partial_derivatives_of_u;
}

template <typename ResultTags, typename DerivativeTags, size_t Dim,
          typename DerivativeFrame>
void partial_derivatives(
    const gsl::not_null<Variables<ResultTags>*> du,
    const std::array<Variables<DerivativeTags>, Dim>&
        logical_partial_derivatives_of_u,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>& inverse_jacobian) {
  using ValueType = typename Variables<DerivativeTags>::value_type;
  auto& partial_derivatives_of_u = *du;
  // For mutating compute items we must set the size.
  if (UNLIKELY(partial_derivatives_of_u.number_of_grid_points() !=
               logical_partial_derivatives_of_u[0].number_of_grid_points())) {
    partial_derivatives_of_u.initialize(
        logical_partial_derivatives_of_u[0].number_of_grid_points());
  }

  std::array<const ValueType*, Dim> logical_derivs{};
  for (size_t i = 0; i < Dim; ++i) {
    gsl::at(logical_derivs, i) =
        gsl::at(logical_partial_derivatives_of_u, i).data();
  }
  partial_derivatives_detail::partial_derivatives_impl(
      make_not_null(&partial_derivatives_of_u), logical_derivs,
      Variables<DerivativeTags>::number_of_independent_components,
      inverse_jacobian);
}

template <typename ResultTags, typename VariableTags, size_t Dim,
          typename DerivativeFrame>
void partial_derivatives(
    const gsl::not_null<Variables<ResultTags>*> du,
    const Variables<VariableTags>& u, const Mesh<Dim>& mesh,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>& inverse_jacobian) {
  using DerivativeTags =
      tmpl::front<tmpl::split_at<VariableTags, tmpl::size<ResultTags>>>;
  static_assert(
      std::is_same_v<
          tmpl::transform<ResultTags, tmpl::bind<tmpl::type_from, tmpl::_1>>,
          tmpl::transform<db::wrap_tags_in<Tags::deriv, DerivativeTags,
                                           tmpl::size_t<Dim>, DerivativeFrame>,
                          tmpl::bind<tmpl::type_from, tmpl::_1>>>);
  using ValueType = typename Variables<VariableTags>::value_type;
  auto& partial_derivatives_of_u = *du;
  // For mutating compute items we must set the size.
  if (UNLIKELY(partial_derivatives_of_u.number_of_grid_points() !=
               mesh.number_of_grid_points())) {
    partial_derivatives_of_u.initialize(mesh.number_of_grid_points());
  }

  const size_t vars_size =
      u.number_of_grid_points() *
      Variables<DerivativeTags>::number_of_independent_components;
  const auto logical_derivs_data =
      cpp20::make_unique_for_overwrite<ValueType[]>(
          (Dim > 1 ? (Dim + 1) : Dim) * vars_size);
  std::array<ValueType*, Dim> logical_derivs{};
  for (size_t i = 0; i < Dim; ++i) {
    gsl::at(logical_derivs, i) = &(logical_derivs_data[i * vars_size]);
  }
  Variables<DerivativeTags> temp{};
  if constexpr (Dim > 1) {
    temp.set_data_ref(&logical_derivs_data[Dim * vars_size], vars_size);
  }
  partial_derivatives_detail::LogicalImpl<
      Dim, VariableTags, DerivativeTags>::apply(make_not_null(&logical_derivs),
                                                &partial_derivatives_of_u,
                                                &temp, u, mesh);

  std::array<const ValueType*, Dim> const_logical_derivs{};
  for (size_t i = 0; i < Dim; ++i) {
    gsl::at(const_logical_derivs, i) = gsl::at(logical_derivs, i);
  }
  partial_derivatives_detail::partial_derivatives_impl(
      make_not_null(&partial_derivatives_of_u), const_logical_derivs,
      Variables<DerivativeTags>::number_of_independent_components,
      inverse_jacobian);
}

template <typename DerivativeTags, typename VariableTags, size_t Dim,
          typename DerivativeFrame>
Variables<db::wrap_tags_in<Tags::deriv, DerivativeTags, tmpl::size_t<Dim>,
                           DerivativeFrame>>
partial_derivatives(
    const Variables<VariableTags>& u, const Mesh<Dim>& mesh,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>& inverse_jacobian) {
  Variables<db::wrap_tags_in<Tags::deriv, DerivativeTags, tmpl::size_t<Dim>,
                             DerivativeFrame>>
      partial_derivatives_of_u(u.number_of_grid_points());
  partial_derivatives(make_not_null(&partial_derivatives_of_u),
                                      u, mesh, inverse_jacobian);
  return partial_derivatives_of_u;
}

namespace partial_derivatives_detail {
template <typename VariableTags, typename DerivativeTags>
struct LogicalImpl<1, VariableTags, DerivativeTags> {
  static constexpr const size_t Dim = 1;
  template <typename T, typename ValueType = typename Variables<T>::value_type>
  static void apply(
      const gsl::not_null<std::array<ValueType*, Dim>*> logical_du,
      Variables<T>* /*unused_in_1d*/,
      Variables<DerivativeTags>* const /*unused_in_1d*/,
      const Variables<VariableTags>& u, const Mesh<Dim>& mesh) {
    auto& logical_partial_derivatives_of_u = *logical_du;
    const size_t deriv_size =
        Variables<DerivativeTags>::number_of_independent_components *
        u.number_of_grid_points();
    const Matrix& differentiation_matrix_xi =
        Spectral::differentiation_matrix(mesh.slice_through(0));
    apply_matrix_in_first_dim(logical_partial_derivatives_of_u[0], u.data(),
                              differentiation_matrix_xi, deriv_size);
  }
};

template <typename VariableTags, typename DerivativeTags>
struct LogicalImpl<2, VariableTags, DerivativeTags> {
  static constexpr size_t Dim = 2;
  template <typename T, typename ValueType = typename Variables<T>::value_type>
  static void apply(
      const gsl::not_null<std::array<ValueType*, Dim>*> logical_du,
      Variables<T>* const partial_u_wrt_eta,
      Variables<DerivativeTags>* const u_eta_fastest,
      const Variables<VariableTags>& u, const Mesh<2>& mesh) {
    static_assert(
        Variables<DerivativeTags>::number_of_independent_components <=
            Variables<T>::number_of_independent_components,
        "Temporary buffer in logical partial derivatives is too small");
    auto& logical_partial_derivatives_of_u = *logical_du;
    const size_t deriv_size =
        Variables<DerivativeTags>::number_of_independent_components *
        u.number_of_grid_points();
    const Matrix& differentiation_matrix_xi =
        Spectral::differentiation_matrix(mesh.slice_through(0));
    const size_t num_components_times_xi_slices = deriv_size / mesh.extents(0);
    apply_matrix_in_first_dim(logical_partial_derivatives_of_u[0], u.data(),
                              differentiation_matrix_xi, deriv_size);
    transpose<Variables<VariableTags>, Variables<DerivativeTags>>(
        make_not_null(u_eta_fastest), u, mesh.extents(0),
        num_components_times_xi_slices);
    const Matrix& differentiation_matrix_eta =
        Spectral::differentiation_matrix(mesh.slice_through(1));
    apply_matrix_in_first_dim(partial_u_wrt_eta->data(), u_eta_fastest->data(),
                              differentiation_matrix_eta, deriv_size);
    raw_transpose(make_not_null(logical_partial_derivatives_of_u[1]),
                  partial_u_wrt_eta->data(), num_components_times_xi_slices,
                  mesh.extents(0));
  }
};

template <typename VariableTags, typename DerivativeTags>
struct LogicalImpl<3, VariableTags, DerivativeTags> {
  static constexpr size_t Dim = 3;
  template <class T, typename ValueType = typename Variables<T>::value_type>
  static void apply(
      const gsl::not_null<std::array<ValueType*, Dim>*> logical_du,
      Variables<T>* const partial_u_wrt_eta_or_zeta,
      Variables<DerivativeTags>* const u_eta_or_zeta_fastest,
      const Variables<VariableTags>& u, const Mesh<3>& mesh) {
    if (mesh.basis(1) == Spectral::Basis::SphericalHarmonic) {
      if constexpr (std::is_same_v<ValueType, double>) {
        spherical_apply(logical_du, u, mesh);
      } else {
        ERROR(
            "Support for complex numbers with spherical harmonics is not yet "
            "implemented for logical_partial_derivative.");
      }
    } else {
      static_assert(
          Variables<DerivativeTags>::number_of_independent_components <=
              Variables<T>::number_of_independent_components,
          "Temporary buffer in logical partial derivatives is too small");
      auto& logical_partial_derivatives_of_u = *logical_du;
      const Matrix& differentiation_matrix_xi =
          Spectral::differentiation_matrix(mesh.slice_through(0));
      const size_t deriv_size =
          Variables<DerivativeTags>::number_of_independent_components *
          u.number_of_grid_points();
      const size_t num_components_times_xi_slices =
          deriv_size / mesh.extents(0);
      apply_matrix_in_first_dim(logical_partial_derivatives_of_u[0], u.data(),
                                differentiation_matrix_xi, deriv_size);
      transpose<Variables<VariableTags>, Variables<DerivativeTags>>(
          make_not_null(u_eta_or_zeta_fastest), u, mesh.extents(0),
          num_components_times_xi_slices);
      const Matrix& differentiation_matrix_eta =
          Spectral::differentiation_matrix(mesh.slice_through(1));
      apply_matrix_in_first_dim(partial_u_wrt_eta_or_zeta->data(),
                                u_eta_or_zeta_fastest->data(),
                                differentiation_matrix_eta, deriv_size);
      raw_transpose(make_not_null(logical_partial_derivatives_of_u[1]),
                    partial_u_wrt_eta_or_zeta->data(),
                    num_components_times_xi_slices, mesh.extents(0));

      const size_t chunk_size = mesh.extents(0) * mesh.extents(1);
      const size_t number_of_chunks = deriv_size / chunk_size;
      transpose(make_not_null(u_eta_or_zeta_fastest), u, chunk_size,
                number_of_chunks);
      const Matrix& differentiation_matrix_zeta =
          Spectral::differentiation_matrix(mesh.slice_through(2));
      apply_matrix_in_first_dim(partial_u_wrt_eta_or_zeta->data(),
                                u_eta_or_zeta_fastest->data(),
                                differentiation_matrix_zeta, deriv_size);
      raw_transpose(make_not_null(logical_partial_derivatives_of_u[2]),
                    partial_u_wrt_eta_or_zeta->data(), number_of_chunks,
                    chunk_size);
    }
  }

  static void spherical_apply(
      const gsl::not_null<std::array<double*, Dim>*> logical_du,
      const Variables<VariableTags>& u, const Mesh<3>& mesh) {
    auto& logical_partial_derivatives_of_u = *logical_du;
    const Matrix& differentiation_matrix_xi =
        Spectral::differentiation_matrix(mesh.slice_through(0));
    const size_t deriv_size =
        Variables<DerivativeTags>::number_of_independent_components *
        u.number_of_grid_points();
    apply_matrix_in_first_dim(logical_partial_derivatives_of_u[0], u.data(),
                              differentiation_matrix_xi, deriv_size);
    const auto& ylm = ylm::get_spherepack_cache(mesh.extents(1) - 1);
    for (size_t n = 0;
         n < Variables<DerivativeTags>::number_of_independent_components; ++n) {
      const size_t offset = n * u.number_of_grid_points();
      const auto du = std::array{logical_partial_derivatives_of_u[1] + offset,
                                 logical_partial_derivatives_of_u[2] + offset};
      ylm.gradient_all_offsets(du, make_not_null(u.data() + offset),
                               mesh.extents(0));
    }
  }
};
}  // namespace partial_derivatives_detail

// Macro to explicitly instantiate partial_derivatives()
// for a given system of equations
#define INSTANTIATE_PARTIAL_DERIVATIVES_WITH_SYSTEM(SYSTEM, DIM,             \
                                                    DERIVATIVE_FRAME)        \
  template void partial_derivatives(                                         \
      gsl::not_null<Variables<                                               \
          db::wrap_tags_in<Tags::deriv, typename SYSTEM::gradient_variables, \
                           tmpl::size_t<DIM>, DERIVATIVE_FRAME>>*>           \
          du,                                                                \
      const std::array<Variables<typename SYSTEM::gradient_variables>, DIM>& \
          logical_partial_derivatives_of_u,                                  \
      const InverseJacobian<DataVector, DIM, Frame::ElementLogical,          \
                            DERIVATIVE_FRAME>& inverse_jacobian);            \
  template void partial_derivatives(                                         \
      gsl::not_null<Variables<                                               \
          db::wrap_tags_in<Tags::deriv, typename SYSTEM::gradient_variables, \
                           tmpl::size_t<DIM>, DERIVATIVE_FRAME>>*>           \
          du,                                                                \
      const Variables<typename SYSTEM::gradient_variables>& u,               \
      const Mesh<DIM>& mesh,                                                 \
      const InverseJacobian<DataVector, DIM, Frame::ElementLogical,          \
                            DERIVATIVE_FRAME>& inverse_jacobian);            \
  template Variables<                                                        \
      db::wrap_tags_in<Tags::deriv, typename SYSTEM::gradient_variables,     \
                       tmpl::size_t<DIM>, DERIVATIVE_FRAME>>                 \
  partial_derivatives<typename SYSTEM::gradient_variables>(                  \
      const Variables<typename SYSTEM::gradient_variables>& u,               \
      const Mesh<DIM>& mesh,                                                 \
      const InverseJacobian<DataVector, DIM, Frame::ElementLogical,          \
                            DERIVATIVE_FRAME>& inverse_jacobian);
