// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include "NumericalAlgorithms/LinearOperators/Divergence.hpp"

#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.tpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Utilities/StdArrayHelpers.hpp"

template <typename FluxTags, size_t Dim, typename DerivativeFrame>
Variables<db::wrap_tags_in<Tags::div, FluxTags>> divergence(
    const Variables<FluxTags>& F, const Mesh<Dim>& mesh,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>& inverse_jacobian) {
  Variables<db::wrap_tags_in<Tags::div, FluxTags>> divergence_of_F(
      F.number_of_grid_points());
  divergence(make_not_null(&divergence_of_F), F, mesh, inverse_jacobian);
  return divergence_of_F;
}

template <typename... DivTags, typename... FluxTags, size_t Dim,
          typename DerivativeFrame>
void divergence(
    const gsl::not_null<Variables<tmpl::list<DivTags...>>*> divergence_of_F,
    const Variables<tmpl::list<FluxTags...>>& F, const Mesh<Dim>& mesh,
    const InverseJacobian<DataVector, Dim, Frame::ElementLogical,
                          DerivativeFrame>& inverse_jacobian) {
  if (UNLIKELY(divergence_of_F->number_of_grid_points() !=
               mesh.number_of_grid_points())) {
    divergence_of_F->initialize(mesh.number_of_grid_points());
  }

  using DerivativeTags = tmpl::list<FluxTags...>;
  using ValueType = typename Variables<DerivativeTags>::value_type;
  const size_t vars_size =
      Variables<DerivativeTags>::number_of_independent_components *
      F.number_of_grid_points();
  const auto logical_derivs_data =
      cpp20::make_unique_for_overwrite<ValueType[]>(
          (Dim > 1 ? (Dim + 2) : Dim) * vars_size);
  std::array<ValueType*, Dim> logical_derivs{};
  std::array<Variables<DerivativeTags>, Dim> logical_partial_derivatives_of_F{};
  for (size_t i = 0; i < Dim; ++i) {
    gsl::at(logical_derivs, i) = &(logical_derivs_data[i * vars_size]);
    gsl::at(logical_partial_derivatives_of_F, i)
        .set_data_ref(gsl::at(logical_derivs, i), vars_size);
  }
  if constexpr (Dim > 1) {
    Variables<DerivativeTags> temp0{&logical_derivs_data[Dim * vars_size],
                                    vars_size};
    Variables<DerivativeTags> temp1{&logical_derivs_data[(Dim + 1) * vars_size],
                                    vars_size};
    partial_derivatives_detail::
        LogicalImpl<Dim, tmpl::list<FluxTags...>, DerivativeTags>::apply(
            make_not_null(&logical_derivs), &temp0, &temp1, F, mesh);
  } else {
    Variables<DerivativeTags> *temp = nullptr;
    partial_derivatives_detail::
        LogicalImpl<Dim, tmpl::list<FluxTags...>, DerivativeTags>::apply(
            make_not_null(&logical_derivs), temp, temp, F, mesh);
  }

  const auto apply_div = [
    &divergence_of_F, &inverse_jacobian, &logical_partial_derivatives_of_F
  ](auto flux_tag_v, auto div_tag_v) {
    using FluxTag = std::decay_t<decltype(flux_tag_v)>;
    using DivFluxTag = std::decay_t<decltype(div_tag_v)>;

    using first_index =
        tmpl::front<typename FluxTag::type::index_list>;
    static_assert(
        std::is_same_v<typename first_index::Frame, DerivativeFrame> and
            first_index::ul == UpLo::Up,
        "First index of tensor cannot be contracted with derivative "
        "because either it is in the wrong frame or it has the wrong "
        "valence");

    auto& divergence_of_flux = get<DivFluxTag>(*divergence_of_F);
    for (auto it = divergence_of_flux.begin(); it != divergence_of_flux.end();
         ++it) {
      *it = 0.0;
      const auto div_flux_indices = divergence_of_flux.get_tensor_index(it);
      for (size_t i0 = 0; i0 < Dim; ++i0) {
        const auto flux_indices = prepend(div_flux_indices, i0);
        for (size_t d = 0; d < Dim; ++d) {
          *it += inverse_jacobian.get(d, i0) *
                 get<FluxTag>(gsl::at(logical_partial_derivatives_of_F, d))
                     .get(flux_indices);
        }
      }
    }
  };
  EXPAND_PACK_LEFT_TO_RIGHT(apply_div(FluxTags{}, DivTags{}));
}

template <typename FluxTags, size_t Dim>
Variables<db::wrap_tags_in<Tags::div, FluxTags>> logical_divergence(
    const Variables<FluxTags>& flux, const Mesh<Dim>& mesh) {
  Variables<db::wrap_tags_in<Tags::div, FluxTags>> div_flux(
      flux.number_of_grid_points());
  logical_divergence(make_not_null(&div_flux), flux, mesh);
  return div_flux;
}

template <typename... ResultTags, typename... FluxTags, size_t Dim>
void logical_divergence(
    const gsl::not_null<Variables<tmpl::list<ResultTags...>>*> div_flux,
    const Variables<tmpl::list<FluxTags...>>& flux, const Mesh<Dim>& mesh) {
  static_assert(
      (... and
       std::is_same_v<tmpl::front<typename FluxTags::type::index_list>,
                      SpatialIndex<Dim, UpLo::Up, Frame::ElementLogical>>),
      "The first index of each flux must be an upper spatial index in "
      "element-logical coordinates.");
  static_assert(
      (... and
       std::is_same_v<
           typename ResultTags::type,
           TensorMetafunctions::remove_first_index<typename FluxTags::type>>),
      "The result tensors must have the same type as the flux tensors with "
      "their first index removed.");
  div_flux->initialize(mesh.number_of_grid_points(), 0.);
  // Note: This function hasn't been optimized much at all. Feel free to
  // optimize if needed!
  EXPAND_PACK_LEFT_TO_RIGHT(logical_divergence(
      make_not_null(&get<ResultTags>(*div_flux)), get<FluxTags>(flux), mesh));
}
