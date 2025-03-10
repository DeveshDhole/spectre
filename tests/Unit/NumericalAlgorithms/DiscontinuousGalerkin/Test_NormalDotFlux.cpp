// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <string>

#include "DataStructures/ComplexDataVector.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Domain/FaceNormal.hpp"
#include "Framework/CheckWithRandomValues.hpp"
#include "Framework/SetupLocalPythonEnvironment.hpp"
#include "Helpers/DataStructures/DataBox/TestHelpers.hpp"
#include "NumericalAlgorithms/DiscontinuousGalerkin/NormalDotFlux.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

namespace {

template <size_t Dim, typename FluxTensor, typename ResultTensor>
void check_normal_dot_flux(const tnsr::i<DataVector, Dim>& normal,
                           const FluxTensor& flux_tensor,
                           const ResultTensor& expected_result) {
  ResultTensor result;
  normal_dot_flux(make_not_null(&result), normal, flux_tensor);
  for (auto it = result.begin(); it != result.end(); it++) {
    CHECK_ITERABLE_APPROX(*it,
                          expected_result.get(result.get_tensor_index(it)));
  }
}

template <size_t Dim, typename Fr, typename DataType, typename Symm,
          typename... RemainingIndices>
void test_with_random_values(
    const DataVector& used_for_size,
    Tensor<DataType, Symm,
           index_list<SpatialIndex<Dim, UpLo::Up, Fr>, RemainingIndices...>>
    /*meta*/) {
  using NDotFluxTensor =
      Tensor<DataType, tmpl::pop_front<Symm>, index_list<RemainingIndices...>>;
  using FluxTensor =
      Tensor<DataType, Symm,
             index_list<SpatialIndex<Dim, UpLo::Up, Fr>, RemainingIndices...>>;
  pypp::check_with_random_values<1>(
      // This static_cast helps GCC figure out the type of the function
      static_cast<void (*)(const gsl::not_null<NDotFluxTensor*>,
                           const tnsr::i<DataVector, Dim, Fr>&,
                           const FluxTensor&)>(
          &normal_dot_flux<Dim, Fr, DataType, Symm, RemainingIndices...>),
      "NormalDotFlux", {"normal_dot_flux"}, {{{-1.0, 1.0}}}, used_for_size);
  pypp::check_with_random_values<1>(
      static_cast<void (*)(
          const gsl::not_null<TensorMetafunctions::prepend_spatial_index<
              NDotFluxTensor, Dim, UpLo::Lo, Fr>*>,
          const tnsr::i<DataVector, Dim, Fr>&, const NDotFluxTensor&)>(
          &normal_times_flux<Dim, Fr, DataType, tmpl::pop_front<Symm>,
                             index_list<RemainingIndices...>>),
      "NormalDotFlux", {"normal_times_flux"}, {{{-1.0, 1.0}}}, used_for_size);
}

struct Var1 : db::SimpleTag {
  using type = Scalar<DataVector>;
};

template <size_t Dim, typename Frame>
struct Var2 : db::SimpleTag {
  using type = tnsr::i<DataVector, Dim, Frame>;
};

template <size_t Dim, typename Frame>
using variables_tag = Tags::Variables<tmpl::list<Var1, Var2<Dim, Frame>>>;

template <size_t Dim, typename Frame>
using flux1 = Tags::Flux<Var1, tmpl::size_t<Dim>, Frame>;
template <size_t Dim, typename Frame>
using flux2 = Tags::Flux<Var2<Dim, Frame>, tmpl::size_t<Dim>, Frame>;

template <size_t Dim, typename Frame>
using flux_tag = db::add_tag_prefix<Tags::Flux, variables_tag<Dim, Frame>,
                                    tmpl::size_t<Dim>, Frame>;

// Copy the values in the components of Tensor<double, ...> `in` into
// the `data_index` entry of DataVectors in the Tensor<DataVector,
// ...> `out`, mapping the component (in0, in1, ...) to (in0, in1,
// ..., `extra_indices`).
template <typename OutTensor, typename InTensor>
void copy_into(const gsl::not_null<OutTensor*> out, const InTensor& in,
               const std::array<size_t, OutTensor::rank() - InTensor::rank()>&
                   extra_indices,
               const size_t data_index) {
  for (auto it = in.begin(); it != in.end(); ++it) {
    const auto in_index = in.get_tensor_index(it);
    std::array<size_t, OutTensor::rank()> out_index{};
    for (size_t i = 0; i < InTensor::rank(); ++i) {
      gsl::at(out_index, i) = gsl::at(in_index, i);
    }
    for (size_t i = 0; i < OutTensor::rank() - InTensor::rank(); ++i) {
      gsl::at(out_index, i + InTensor::rank()) = gsl::at(extra_indices, i);
    }
    out->get(out_index)[data_index] = in.get(in_index);
  }
}

template <size_t Dim, typename Frame>
tnsr::i<double, Dim, Frame> generate_normal(const size_t seed) {
  tnsr::i<double, Dim, Frame> result{};
  std::iota(result.begin(), result.end(), seed + 2.);
  return result;
}

template <size_t Dim, typename Frame>
tnsr::I<double, Dim, Frame> generate_flux(const size_t seed) {
  tnsr::I<double, Dim, Frame> result{};
  std::iota(result.begin(), result.end(), seed + 3.);
  return result;
}

template <size_t Dim>
Scalar<double> generate_f_dot_n(const size_t normal_seed,
                                const size_t flux_seed) {
  double magnitude_normal = 0.;
  double unnormalized_f_dot_n = 0.;
  for (size_t i = 0; i < Dim; ++i) {
    magnitude_normal += square(static_cast<double>(normal_seed + i) + 2.);
    unnormalized_f_dot_n += (static_cast<double>(normal_seed + i) + 2.) *
                            (static_cast<double>(flux_seed + i) + 3.);
  }
  magnitude_normal = sqrt(magnitude_normal);

  return Scalar<double>(unnormalized_f_dot_n / magnitude_normal);
}

template <size_t Dim>
void test_with_variables() {
  using Fr = Frame::Inertial;
  constexpr size_t num_points = 5;
  tnsr::i<DataVector, Dim, Fr> normal(num_points);
  typename flux_tag<Dim, Fr>::type fluxes(num_points);
  Var1::type expected1(num_points);
  typename Var2<Dim, Fr>::type expected2(num_points);
  for (size_t i = 0; i < num_points; ++i) {
    copy_into(make_not_null(&normal), generate_normal<Dim, Fr>(i), {}, i);
    copy_into(make_not_null(&get<flux1<Dim, Fr>>(fluxes)),
              generate_flux<Dim, Fr>(i), {}, i);
    copy_into(make_not_null(&expected1), generate_f_dot_n<Dim>(i, i), {}, i);
    for (size_t j = 0; j < Dim; ++j) {
      copy_into(make_not_null(&get<flux2<Dim, Fr>>(fluxes)),
                generate_flux<Dim, Fr>(i + 10 * j), {{j}}, i);
      copy_into(make_not_null(&expected2), generate_f_dot_n<Dim>(i, i + 10 * j),
                {{j}}, i);
    }
  }

  const auto magnitude_normal = magnitude(normal);
  for (size_t d = 0; d < Dim; d++) {
    normal.get(d) /= get(magnitude_normal);
  }

  const auto result =
      normal_dot_flux<typename variables_tag<Dim, Fr>::tags_list>(normal,
                                                                  fluxes);

  CHECK_ITERABLE_APPROX(get<Tags::NormalDotFlux<Var1>>(result), expected1);
  CHECK_ITERABLE_APPROX((get<Tags::NormalDotFlux<Var2<Dim, Fr>>>(result)),
                        expected2);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.DiscontinuousGalerkin.NormalDotFlux",
                  "[Unit][Evolution]") {
  {
    INFO("Explicit values");
    const size_t npts = 5;
    const DataVector zero(npts, 0.0);
    const DataVector one(npts, 1.0);
    const DataVector two(npts, 2.0);
    const DataVector three(npts, 3.0);
    const DataVector four(npts, 4.0);
    const DataVector five(npts, 5.0);
    const DataVector six(npts, 6.0);
    const DataVector seven(npts, 7.0);
    const DataVector eight(npts, 8.0);
    const DataVector nine(npts, 9.0);
    const DataVector ten(npts, 10.0);
    const DataVector eleven(npts, 11.0);
    const DataVector twelve(npts, 12.0);
    const DataVector fifteen(npts, 15.0);
    const DataVector eighteen(npts, 18.0);

    check_normal_dot_flux(tnsr::i<DataVector, 1>{{{one}}},
                          tnsr::I<DataVector, 1>{{{two}}},
                          Scalar<DataVector>{two});
    check_normal_dot_flux(tnsr::i<DataVector, 2>{{{one, two}}},
                          tnsr::I<DataVector, 2>{{{three, four}}},
                          Scalar<DataVector>{eleven});
    check_normal_dot_flux(tnsr::i<DataVector, 3>{{{one, two, three}}},
                          tnsr::I<DataVector, 3>{{{-four, -two, three}}},
                          Scalar<DataVector>{one});
    check_normal_dot_flux(tnsr::i<DataVector, 2>{{{one, two}}},
                          [&one, &two, &three, &four]() {
                            tnsr::Ij<DataVector, 2> flux;
                            get<0, 0>(flux) = one;
                            get<0, 1>(flux) = two;
                            get<1, 0>(flux) = three;
                            get<1, 1>(flux) = four;
                            return flux;
                          }(),
                          tnsr::i<DataVector, 2>{{{seven, ten}}});
    check_normal_dot_flux(tnsr::i<DataVector, 2>{{{one, two}}},
                          [&one, &two, &three]() {
                            tnsr::II<DataVector, 2> flux;
                            get<0, 0>(flux) = one;
                            get<0, 1>(flux) = two;
                            get<1, 1>(flux) = three;
                            return flux;
                          }(),
                          tnsr::I<DataVector, 2>{{{five, eight}}});
    check_normal_dot_flux(tnsr::i<DataVector, 3>{{{one, two, three}}},
                          [&one, &two, &three, &four, &five, &six]() {
                            tnsr::II<DataVector, 3> flux;
                            get<0, 0>(flux) = one;
                            get<0, 1>(flux) = two;
                            get<0, 2>(flux) = -three;
                            get<1, 1>(flux) = -four;
                            get<1, 2>(flux) = five;
                            get<2, 2>(flux) = -six;
                            return flux;
                          }(),
                          tnsr::I<DataVector, 3>{{{-four, nine, -eleven}}});
    check_normal_dot_flux(tnsr::i<DataVector, 2>{{{one, two}}},
                          [&one, &two, &three, &four, &five, &six]() {
                            tnsr::Iaa<DataVector, 2> flux;
                            get<0, 0, 0>(flux) = one;
                            get<0, 0, 1>(flux) = two;
                            get<0, 0, 2>(flux) = -three;
                            get<0, 1, 1>(flux) = -four;
                            get<0, 1, 2>(flux) = five;
                            get<0, 2, 2>(flux) = -six;
                            get<1, 0, 0>(flux) = two;
                            get<1, 0, 1>(flux) = one;
                            get<1, 0, 2>(flux) = -three;
                            get<1, 1, 1>(flux) = two;
                            get<1, 1, 2>(flux) = three;
                            get<1, 2, 2>(flux) = three;
                            return flux;
                          }(),
                          [&five, &four, &nine, &zero, &eleven]() {
                            tnsr::aa<DataVector, 2> result;
                            get<0, 0>(result) = five;
                            get<0, 1>(result) = four;
                            get<0, 2>(result) = -nine;
                            get<1, 1>(result) = zero;
                            get<1, 2>(result) = eleven;
                            get<2, 2>(result) = zero;
                            return result;
                          }());
    check_normal_dot_flux(
        tnsr::i<DataVector, 2>{{{one, two}}},
        [&one, &two, &three, &four, &five, &six]() {
          tnsr::Ijaa<DataVector, 2> flux;
          get<0, 0, 0, 0>(flux) = one;
          get<0, 0, 0, 1>(flux) = two;
          get<0, 0, 0, 2>(flux) = -three;
          get<0, 0, 1, 1>(flux) = -four;
          get<0, 0, 1, 2>(flux) = five;
          get<0, 0, 2, 2>(flux) = -six;
          get<0, 1, 0, 0>(flux) = two;
          get<0, 1, 0, 1>(flux) = one;
          get<0, 1, 0, 2>(flux) = -three;
          get<0, 1, 1, 1>(flux) = two;
          get<0, 1, 1, 2>(flux) = three;
          get<0, 1, 2, 2>(flux) = three;
          get<1, 0, 0, 0>(flux) = one;
          get<1, 0, 0, 1>(flux) = two;
          get<1, 0, 0, 2>(flux) = -three;
          get<1, 0, 1, 1>(flux) = -four;
          get<1, 0, 1, 2>(flux) = five;
          get<1, 0, 2, 2>(flux) = -six;
          get<1, 1, 0, 0>(flux) = two;
          get<1, 1, 0, 1>(flux) = one;
          get<1, 1, 0, 2>(flux) = -three;
          get<1, 1, 1, 1>(flux) = two;
          get<1, 1, 1, 2>(flux) = three;
          get<1, 1, 2, 2>(flux) = three;
          return flux;
        }(),
        [&three, &six, &nine, &twelve, &fifteen, &eighteen]() {
          tnsr::iaa<DataVector, 2> result;
          get<0, 0, 0>(result) = three;
          get<0, 0, 1>(result) = six;
          get<0, 0, 2>(result) = -nine;
          get<0, 1, 1>(result) = -twelve;
          get<0, 1, 2>(result) = fifteen;
          get<0, 2, 2>(result) = -eighteen;
          get<1, 0, 0>(result) = six;
          get<1, 0, 1>(result) = three;
          get<1, 0, 2>(result) = -nine;
          get<1, 1, 1>(result) = six;
          get<1, 1, 2>(result) = nine;
          get<1, 2, 2>(result) = nine;
          return result;
        }());
  }
  {
    INFO("Random values");
    pypp::SetupLocalPythonEnvironment local_python_env{"Domain"};

    GENERATE_UNINITIALIZED_DATAVECTOR;
    test_with_random_values(dv, tnsr::I<DataVector, 1>{});
    test_with_random_values(dv, tnsr::I<DataVector, 2>{});
    test_with_random_values(dv, tnsr::I<DataVector, 3>{});
    test_with_random_values(dv, tnsr::I<ComplexDataVector, 3>{});
    test_with_random_values(dv, tnsr::II<DataVector, 1>{});
    test_with_random_values(dv, tnsr::II<DataVector, 2>{});
    test_with_random_values(dv, tnsr::II<DataVector, 3>{});
    test_with_random_values(dv, tnsr::Ij<DataVector, 1>{});
    test_with_random_values(dv, tnsr::Ij<DataVector, 2>{});
    test_with_random_values(dv, tnsr::Ij<DataVector, 3>{});
    test_with_random_values(dv, tnsr::Ij<ComplexDataVector, 3>{});
    test_with_random_values(dv, tnsr::Ijk<DataVector, 1>{});
    test_with_random_values(dv, tnsr::Ijk<DataVector, 2>{});
    test_with_random_values(dv, tnsr::Ijk<DataVector, 3>{});
    test_with_random_values(dv, tnsr::III<DataVector, 1>{});
    test_with_random_values(dv, tnsr::III<DataVector, 2>{});
    test_with_random_values(dv, tnsr::III<DataVector, 3>{});
    test_with_random_values(dv, tnsr::Iaa<DataVector, 1>{});
    test_with_random_values(dv, tnsr::Iaa<DataVector, 2>{});
    test_with_random_values(dv, tnsr::Iaa<DataVector, 3>{});
    test_with_random_values(dv, tnsr::Ijaa<DataVector, 1>{});
    test_with_random_values(dv, tnsr::Ijaa<DataVector, 2>{});
    test_with_random_values(dv, tnsr::Ijaa<DataVector, 3>{});
  }
  {
    INFO("Variables");
    GENERATE_UNINITIALIZED_DATAVECTOR;
    test_with_variables<1>();
    test_with_variables<2>();
    test_with_variables<3>();
  }
}
