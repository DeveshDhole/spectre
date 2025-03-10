// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/AnalyticData/GrMhd/SpecInitialData.hpp"

#include <Exporter.hpp>  // The SpEC Exporter
#include <memory>
#include <optional>
#include <pup.h>
#include <string>
#include <utility>
#include <vector>

#include "DataStructures/Tensor/EagerMath/DeterminantAndInverse.hpp"
#include "DataStructures/Tensor/EagerMath/DotProduct.hpp"
#include "DataStructures/Tensor/EagerMath/RaiseOrLowerIndex.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "IO/External/InterpolateFromSpec.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/Hydro/SpecificEnthalpy.hpp"
#include "PointwiseFunctions/Hydro/Tags.hpp"
#include "Utilities/ContainerHelpers.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/GenerateInstantiations.hpp"
#include "Utilities/Serialization/PupStlCpp17.hpp"
#include "Utilities/System/ParallelInfo.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace grmhd::AnalyticData {

template <size_t ThermodynamicDim>
SpecInitialData<ThermodynamicDim>::SpecInitialData(
    std::string data_directory,
    std::unique_ptr<equation_of_state_type> equation_of_state,
    const double density_cutoff, const double atmosphere_density,
    const std::optional<double> electron_fraction)
    : data_directory_(std::move(data_directory)),
      equation_of_state_(std::move(equation_of_state)),
      density_cutoff_(density_cutoff),
      atmosphere_density_(atmosphere_density),
      electron_fraction_(electron_fraction),
      spec_exporter_(std::make_unique<spec::Exporter>(
          sys::procs_on_node(sys::my_node()), data_directory_,
          vars_to_interpolate_)) {}

template <size_t ThermodynamicDim>
SpecInitialData<ThermodynamicDim>::SpecInitialData(const SpecInitialData& rhs)
    : evolution::initial_data::InitialData(rhs) {
  *this = rhs;
}

template <size_t ThermodynamicDim>
SpecInitialData<ThermodynamicDim>& SpecInitialData<ThermodynamicDim>::operator=(
    const SpecInitialData& rhs) {
  data_directory_ = rhs.data_directory_;
  equation_of_state_ = rhs.equation_of_state_->get_clone();
  density_cutoff_ = rhs.density_cutoff_;
  atmosphere_density_ = rhs.atmosphere_density_;
  electron_fraction_ = rhs.electron_fraction_;
  spec_exporter_ =
      std::make_unique<spec::Exporter>(sys::procs_on_node(sys::my_node()),
                                       data_directory_, vars_to_interpolate_);
  return *this;
}

template <size_t ThermodynamicDim>
std::unique_ptr<evolution::initial_data::InitialData>
SpecInitialData<ThermodynamicDim>::get_clone() const {
  return std::make_unique<SpecInitialData>(*this);
}

template <size_t ThermodynamicDim>
SpecInitialData<ThermodynamicDim>::SpecInitialData(CkMigrateMessage* msg)
    : InitialData(msg) {}

template <size_t ThermodynamicDim>
void SpecInitialData<ThermodynamicDim>::pup(PUP::er& p) {
  InitialData::pup(p);
  p | data_directory_;
  p | equation_of_state_;
  p | density_cutoff_;
  p | atmosphere_density_;
  p | electron_fraction_;
  if (p.isUnpacking()) {
    spec_exporter_ =
        std::make_unique<spec::Exporter>(sys::procs_on_node(sys::my_node()),
                                         data_directory_, vars_to_interpolate_);
  }
}

template <size_t ThermodynamicDim>
PUP::able::PUP_ID SpecInitialData<ThermodynamicDim>::my_PUP_ID = 0;

template <size_t ThermodynamicDim>
template <typename DataType>
tuples::tagged_tuple_from_typelist<typename SpecInitialData<
    ThermodynamicDim>::template interpolated_tags<DataType>>
SpecInitialData<ThermodynamicDim>::interpolate_from_spec(
    const tnsr::I<DataType, 3>& x) const {
  auto interpolated_data =
      io::interpolate_from_spec<interpolated_tags<DataType>>(
          make_not_null(spec_exporter_.get()), x,
          static_cast<size_t>(sys::my_local_rank()));
  DataType& rest_mass_density =
      get(get<hydro::Tags::RestMassDensity<DataType>>(interpolated_data));
  for (size_t i = 0; i < get_size(rest_mass_density); ++i) {
    if (get_element(rest_mass_density, i) < density_cutoff_) {
      get_element(rest_mass_density, i) = atmosphere_density_;
    }
  }
  return interpolated_data;
}

template <size_t ThermodynamicDim>
template <typename DataType>
void SpecInitialData<ThermodynamicDim>::VariablesComputer<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> specific_internal_energy,
    const gsl::not_null<Cache*> cache,
    hydro::Tags::SpecificInternalEnergy<DataType> /*meta*/) const {
  const auto& rest_mass_density =
      get<hydro::Tags::RestMassDensity<DataType>>(interpolated_data);
  const auto& temperature =
      cache->get_var(*this, hydro::Tags::Temperature<DataType>{});
  const auto& electron_fraction =
      cache->get_var(*this, hydro::Tags::ElectronFraction<DataType>{});
  const size_t num_points = get_size(get(rest_mass_density));
  for (size_t i = 0; i < num_points; ++i) {
    const double local_rest_mass_density =
        get_element(get(rest_mass_density), i);
    const double local_temperature = get_element(get(temperature), i);
    const double local_electron_fraction =
        get_element(get(electron_fraction), i);
    if constexpr (ThermodynamicDim == 1) {
      get_element(get(*specific_internal_energy), i) =
          get(eos.specific_internal_energy_from_density(
              Scalar<double>(local_rest_mass_density)));
    } else if constexpr (ThermodynamicDim == 2) {
      get_element(get(*specific_internal_energy), i) =
          get(eos.specific_internal_energy_from_density_and_temperature(
              Scalar<double>(local_rest_mass_density),
              Scalar<double>(local_temperature)));
    } else {
      static_assert(ThermodynamicDim == 3);
      get_element(get(*specific_internal_energy), i) =
          get(eos.specific_internal_energy_from_density_and_temperature(
              Scalar<double>(local_rest_mass_density),
              Scalar<double>(local_temperature),
              Scalar<double>(local_electron_fraction)));
    }
  }
}

template <size_t ThermodynamicDim>
template <typename DataType>
void SpecInitialData<ThermodynamicDim>::VariablesComputer<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> pressure,
    const gsl::not_null<Cache*> cache,
    hydro::Tags::Pressure<DataType> /*meta*/) const {
  const auto& rest_mass_density =
      get<hydro::Tags::RestMassDensity<DataType>>(interpolated_data);
  const auto& temperature =
      cache->get_var(*this, hydro::Tags::Temperature<DataType>{});
  const auto& electron_fraction =
      cache->get_var(*this, hydro::Tags::ElectronFraction<DataType>{});
  const size_t num_points = get_size(get(rest_mass_density));
  for (size_t i = 0; i < num_points; ++i) {
    const double local_rest_mass_density =
        get_element(get(rest_mass_density), i);
    const double local_temperature = get_element(get(temperature), i);
    const double local_electron_fraction =
        get_element(get(electron_fraction), i);
    if constexpr (ThermodynamicDim == 1) {
      get_element(get(*pressure), i) = get(
          eos.pressure_from_density(Scalar<double>(local_rest_mass_density)));
    } else if constexpr (ThermodynamicDim == 2) {
      get_element(get(*pressure), i) = get(eos.pressure_from_density_and_energy(
          Scalar<double>(local_rest_mass_density),
          Scalar<double>(
              get(eos.specific_internal_energy_from_density_and_temperature(
                  Scalar<double>(local_rest_mass_density),
                  Scalar<double>(local_temperature))))));
    } else {
      static_assert(ThermodynamicDim == 3);
      get_element(get(*pressure), i) =
          get(eos.pressure_from_density_and_temperature(
              Scalar<double>(local_rest_mass_density),
              Scalar<double>(local_temperature),
              Scalar<double>(local_electron_fraction)));
    }
  }
}

template <size_t ThermodynamicDim>
template <typename DataType>
void SpecInitialData<ThermodynamicDim>::VariablesComputer<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> specific_enthalpy,
    const gsl::not_null<Cache*> cache,
    hydro::Tags::SpecificEnthalpy<DataType> /*meta*/) const {
  const auto& rest_mass_density =
      get<hydro::Tags::RestMassDensity<DataType>>(interpolated_data);
  const auto& pressure =
      cache->get_var(*this, hydro::Tags::Pressure<DataType>{});
  const auto& specific_internal_energy =
      cache->get_var(*this, hydro::Tags::SpecificInternalEnergy<DataType>{});
  const size_t num_points = get_size(get(rest_mass_density));
  for (size_t i = 0; i < num_points; ++i) {
    const double local_rest_mass_density =
        get_element(get(rest_mass_density), i);
    if (local_rest_mass_density <= density_cutoff) {
      get_element(get(*specific_enthalpy), i) = 1.;
    } else {
      get_element(get(*specific_enthalpy), i) =
          get(hydro::relativistic_specific_enthalpy(
              Scalar<double>(local_rest_mass_density),
              Scalar<double>(get_element(get(specific_internal_energy), i)),
              Scalar<double>(get_element(get(pressure), i))));
    }
  }
}

template <size_t ThermodynamicDim>
template <typename DataType>
void SpecInitialData<ThermodynamicDim>::VariablesComputer<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> temperature,
    const gsl::not_null<Cache*> /*cache*/,
    hydro::Tags::Temperature<DataType> /*meta*/) const {
  const auto& rest_mass_density =
      get<hydro::Tags::RestMassDensity<DataType>>(interpolated_data);
  const size_t num_points = get_size(get(rest_mass_density));
  // Consistent with the use of temperature above we set T=0
  for (size_t i = 0; i < num_points; ++i) {
    get_element(get(*temperature), i) = eos.temperature_lower_bound();
  }
}

template <size_t ThermodynamicDim>
template <typename DataType>
void SpecInitialData<ThermodynamicDim>::VariablesComputer<DataType>::operator()(
    const gsl::not_null<tnsr::II<DataType, 3>*> inv_spatial_metric,
    const gsl::not_null<Cache*> /*cache*/,
    gr::Tags::InverseSpatialMetric<DataType, 3> /*meta*/) const {
  const auto& spatial_metric =
      get<gr::Tags::SpatialMetric<DataType, 3>>(interpolated_data);
  Scalar<DataType> unused_det{};
  determinant_and_inverse(make_not_null(&unused_det), inv_spatial_metric,
                          spatial_metric);
}

template <size_t ThermodynamicDim>
template <typename DataType>
void SpecInitialData<ThermodynamicDim>::VariablesComputer<DataType>::operator()(
    const gsl::not_null<tnsr::I<DataType, 3>*>
        lorentz_factor_times_spatial_velocity,
    const gsl::not_null<Cache*> cache,
    hydro::Tags::LorentzFactorTimesSpatialVelocity<DataType, 3> /*meta*/)
    const {
  const auto& u_i = get<hydro::Tags::LowerSpatialFourVelocity<DataType, 3>>(
      interpolated_data);
  const auto& inv_spatial_metric =
      cache->get_var(*this, gr::Tags::InverseSpatialMetric<DataType, 3>{});
  raise_or_lower_index(lorentz_factor_times_spatial_velocity, u_i,
                       inv_spatial_metric);
}

template <size_t ThermodynamicDim>
template <typename DataType>
void SpecInitialData<ThermodynamicDim>::VariablesComputer<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> lorentz_factor,
    const gsl::not_null<Cache*> cache,
    hydro::Tags::LorentzFactor<DataType> /*meta*/) const {
  const auto& u_i = get<hydro::Tags::LowerSpatialFourVelocity<DataType, 3>>(
      interpolated_data);
  const auto& lorentz_factor_times_spatial_velocity = cache->get_var(
      *this, hydro::Tags::LorentzFactorTimesSpatialVelocity<DataType, 3>{});
  dot_product(lorentz_factor, u_i, lorentz_factor_times_spatial_velocity);
  get(*lorentz_factor) = sqrt(1.0 + get(*lorentz_factor));
}

template <size_t ThermodynamicDim>
template <typename DataType>
void SpecInitialData<ThermodynamicDim>::VariablesComputer<DataType>::operator()(
    const gsl::not_null<tnsr::I<DataType, 3>*> spatial_velocity,
    const gsl::not_null<Cache*> cache,
    hydro::Tags::SpatialVelocity<DataType, 3> /*meta*/) const {
  *spatial_velocity = cache->get_var(
      *this, hydro::Tags::LorentzFactorTimesSpatialVelocity<DataType, 3>{});
  const auto& lorentz_factor =
      cache->get_var(*this, hydro::Tags::LorentzFactor<DataType>{});
  for (size_t d = 0; d < 3; ++d) {
    spatial_velocity->get(d) /= get(lorentz_factor);
  }
}

template <size_t ThermodynamicDim>
template <typename DataType>
void SpecInitialData<ThermodynamicDim>::VariablesComputer<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> electron_fraction,
    const gsl::not_null<Cache*> cache,
    hydro::Tags::ElectronFraction<DataType> /*meta*/) const {
  const auto& temperature =
      cache->get_var(*this, hydro::Tags::Temperature<DataType>{});
  const auto& rest_mass_density =
      get<hydro::Tags::RestMassDensity<DataType>>(interpolated_data);
  const size_t num_points = get_size(get(rest_mass_density));
  if (electron_fraction_value.has_value()) {
    std::fill(electron_fraction->begin(), electron_fraction->end(),
              electron_fraction_value.value());
  } else {
    for (size_t i = 0; i < num_points; ++i) {
      const double local_rest_mass_density =
          get_element(get(rest_mass_density), i);
      const double local_temperature = get_element(get(temperature), i);
      get_element(get(*electron_fraction), i) =
          get(eos.equilibrium_electron_fraction_from_density_temperature(
              Scalar<double>(local_rest_mass_density),
              Scalar<double>(local_temperature)));
    }
  }
}

template <size_t ThermodynamicDim>
template <typename DataType>
void SpecInitialData<ThermodynamicDim>::VariablesComputer<DataType>::operator()(
    const gsl::not_null<tnsr::I<DataType, 3>*> magnetic_field,
    const gsl::not_null<Cache*> /*cache*/,
    hydro::Tags::MagneticField<DataType, 3> /*meta*/) const {
  std::fill(magnetic_field->begin(), magnetic_field->end(), 0.);
}

template <size_t ThermodynamicDim>
template <typename DataType>
void SpecInitialData<ThermodynamicDim>::VariablesComputer<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> div_cleaning_field,
    const gsl::not_null<Cache*> /*cache*/,
    hydro::Tags::DivergenceCleaningField<DataType> /*meta*/) const {
  get(*div_cleaning_field) = 0.;
}

#define THERMODIM(data) BOOST_PP_TUPLE_ELEM(0, data)

#define INSTANTIATION(r, data)                                                \
  template class SpecInitialData<THERMODIM(data)>;                            \
  template tuples::tagged_tuple_from_typelist<typename SpecInitialData<       \
      THERMODIM(data)>::template interpolated_tags<double>>                   \
  SpecInitialData<THERMODIM(data)>::interpolate_from_spec(                    \
      const tnsr::I<double, 3>& x) const;                                     \
  template tuples::tagged_tuple_from_typelist<typename SpecInitialData<       \
      THERMODIM(data)>::template interpolated_tags<DataVector>>               \
  SpecInitialData<THERMODIM(data)>::interpolate_from_spec(                    \
      const tnsr::I<DataVector, 3>& x) const;                                 \
  template class SpecInitialData<THERMODIM(data)>::VariablesComputer<double>; \
  template class SpecInitialData<THERMODIM(                                   \
      data)>::VariablesComputer<DataVector>;

GENERATE_INSTANTIATIONS(INSTANTIATION, (1, 2, 3))

#undef INSTANTIATION
#undef THERMODIM
}  // namespace grmhd::AnalyticData
