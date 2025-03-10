// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <Exporter.hpp>  // The SpEC Exporter
#include <optional>
#include <memory>

#include "DataStructures/CachedTempBuffer.hpp"
#include "Evolution/NumericInitialData.hpp"
#include "Options/String.hpp"
#include "PointwiseFunctions/AnalyticData/GrMhd/AnalyticData.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/Hydro/EquationsOfState/EquationOfState.hpp"
#include "PointwiseFunctions/Hydro/Tags.hpp"
#include "PointwiseFunctions/InitialDataUtilities/InitialData.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

/// \cond
namespace PUP {
class er;
}  // namespace PUP
/// \endcond

namespace grmhd::AnalyticData {

/*!
 * \brief Hydro initial data generated by SpEC.
 *
 * This class loads numerical data written out by the SpEC initial data solver.
 * It uses the `spec::Exporter` linked in from SpEC to interpolate to arbitrary
 * grid points. The coordinates are assumed to be in SpEC's "grid" frame.
 * We interpolate the following quantities:
 *
 * - "g": spatial metric
 * - "K": (lower) extrinsic curvature
 * - "Lapse": lapse
 * - "Shift": (upper) shift
 * - "BaryonDensity": rest mass density
 * - "u_i": lower spatial four-velocity
 *
 * The remaining hydro quantities are computed from the interpolated data and
 * the equation of state.
 * The magnetic field is set to zero and the electron fraction is set to a
 * constant read from the input file.
 */
template <size_t ThermodynamicDim>
class SpecInitialData : public evolution::initial_data::InitialData,
                        public evolution::NumericInitialData,
                        public AnalyticDataBase {
 private:
  struct FromEos {};

 public:
  using equation_of_state_type =
      EquationsOfState::EquationOfState<true, ThermodynamicDim>;

  template <typename DataType>
  using tags = tmpl::append<
      tmpl::list<gr::Tags::SpatialMetric<DataType, 3>,
                 gr::Tags::ExtrinsicCurvature<DataType, 3>,
                 gr::Tags::Lapse<DataType>, gr::Tags::Shift<DataType, 3>>,
      hydro::grmhd_tags<DataType>>;

  static std::string name() {
    return "SpecInitialData" + std::to_string(ThermodynamicDim) + "dEos";
  }

  struct DataDirectory {
    using type = std::string;
    static constexpr Options::String help = {
        "Path to a directory of data produced by SpEC. The directory is "
        "expected to contain 'GrDomain.input' and 'Vars*.h5' files for all the "
        "subdomains in GrDomain.input."};
  };

  struct DensityCutoff {
    using type = double;
    static constexpr Options::String help =
        "Where the density is below this cutoff the fluid variables are set to "
        "vacuum (atmosphere density, atmosphere pressure, atmosphere energy "
        "density/temperature and velocity, unit Lorentz factor and enthalpy). "
        "During the evolution, atmosphere treatment will typically kick in and "
        "fix the value of the fluid variables in these regions. Therefore, "
        "it makes sense to set this density cutoff to the same value as the "
        "atmosphere density cutoff.";
  };

  struct AtmosphereDensity {
    using type = double;
    static constexpr Options::String help =
        "The value to set the density to in atmosphere. This must match what "
        "is done during the evolution.";
  };

  struct ElectronFraction {
    using type = Options::Auto<double, FromEos>;
    static constexpr Options::String help = {
        "Constant electron fraction.\n"
        "To calculate electron fraction from the Eos, set to 'FromEos'."};
  };

  using options = tmpl::list<
      DataDirectory,
      hydro::OptionTags::InitialDataEquationOfState<true, ThermodynamicDim>,
      DensityCutoff, AtmosphereDensity, ElectronFraction>;

  static constexpr Options::String help = {"Initial data generated by SpEC"};

  SpecInitialData() = default;
  SpecInitialData(const SpecInitialData& rhs);
  SpecInitialData& operator=(const SpecInitialData& rhs);
  SpecInitialData(SpecInitialData&& /*rhs*/) = default;
  SpecInitialData& operator=(SpecInitialData&& /*rhs*/) = default;
  ~SpecInitialData() override = default;

  SpecInitialData(std::string data_directory,
                  std::unique_ptr<equation_of_state_type> equation_of_state,
                  double density_cutoff, double atmosphere_density,
                  std::optional<double> electron_fraction);

  auto get_clone() const
      -> std::unique_ptr<evolution::initial_data::InitialData> override;

  /// \cond
  explicit SpecInitialData(CkMigrateMessage* msg);
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(SpecInitialData);
  /// \endcond

  const equation_of_state_type& equation_of_state() const {
    return *equation_of_state_;
  }

  template <typename DataType, typename... Tags>
  tuples::TaggedTuple<Tags...> variables(const tnsr::I<DataType, 3>& x,
                                         tmpl::list<Tags...> /*meta*/) const {
    auto interpolated_vars = interpolate_from_spec(x);
    using requested_tags = tmpl::list<Tags...>;
    using requested_derived_tags =
        tmpl::list_difference<requested_tags, interpolated_tags<DataType>>;
    using requested_interpolated_tags =
        tmpl::list_difference<requested_tags, requested_derived_tags>;
    tuples::TaggedTuple<Tags...> result{};
    // First, compute derived quantities from interpolated data
    if constexpr (tmpl::size<requested_derived_tags>::value > 0) {
      VariablesCache<DataType> cache{get_size(get<0>(x))};
      const VariablesComputer<DataType> computer{
          interpolated_vars, *equation_of_state_, density_cutoff_,
          atmosphere_density_, electron_fraction_};
      tmpl::for_each<requested_derived_tags>(
          [&result, &cache, &computer](const auto tag_v) {
            using tag = tmpl::type_from<std::decay_t<decltype(tag_v)>>;
            get<tag>(result) = cache.get_var(computer, tag{});
          });
    }
    // Then, move interpolated data into result buffer
    tmpl::for_each<requested_interpolated_tags>(
        [&result, &interpolated_vars](const auto tag_v) {
          using tag = tmpl::type_from<std::decay_t<decltype(tag_v)>>;
          get<tag>(result) = std::move(get<tag>(interpolated_vars));
        });
    return result;
  }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& /*p*/) override;

 private:
  /// These quantities are supported for interpolation from SpEC
  template <typename DataType>
  using interpolated_tags = tmpl::list<
      // GR quantities
      gr::Tags::SpatialMetric<DataType, 3>,
      gr::Tags::ExtrinsicCurvature<DataType, 3>, gr::Tags::Lapse<DataType>,
      gr::Tags::Shift<DataType, 3>,
      // Hydro quantities
      hydro::Tags::RestMassDensity<DataType>,
      hydro::Tags::LowerSpatialFourVelocity<DataType, 3>>;

  /// These are the names in SpEC datasets corresponding to the quantities above
  static const inline std::vector<std::string> vars_to_interpolate_{
      // GR quantities
      "g", "K", "Lapse", "Shift",
      // Hydro quantities
      "BaryonDensity", "u_i"};

  template <typename DataType>
  tuples::tagged_tuple_from_typelist<interpolated_tags<DataType>>
  interpolate_from_spec(const tnsr::I<DataType, 3>& x) const;

  /// This cache computes all derived quantities from the interpolated
  /// quantities on demand
  template <typename DataType>
  using VariablesCache = cached_temp_buffer_from_typelist<tmpl::push_back<
      tmpl::list_difference<tags<DataType>, interpolated_tags<DataType>>,
      hydro::Tags::LorentzFactorTimesSpatialVelocity<DataType, 3>,
      gr::Tags::InverseSpatialMetric<DataType, 3>>>;

  /// This is the computer for the `VariablesCache`
  template <typename DataType>
  struct VariablesComputer {
    using Cache = VariablesCache<DataType>;

    const tuples::tagged_tuple_from_typelist<interpolated_tags<DataType>>
        interpolated_data;
    const equation_of_state_type& eos;
    const double density_cutoff;
    const double atmosphere_density;
    const std::optional<double>& electron_fraction_value;

    void operator()(
        gsl::not_null<Scalar<DataType>*> specific_internal_energy,
        gsl::not_null<Cache*> cache,
        hydro::Tags::SpecificInternalEnergy<DataType> /*meta*/) const;
    void operator()(gsl::not_null<Scalar<DataType>*> pressure,
                    gsl::not_null<Cache*> cache,
                    hydro::Tags::Pressure<DataType> /*meta*/) const;
    void operator()(gsl::not_null<Scalar<DataType>*> specific_enthalpy,
                    gsl::not_null<Cache*> cache,
                    hydro::Tags::SpecificEnthalpy<DataType> /*meta*/) const;
    void operator()(gsl::not_null<Scalar<DataType>*> temperature,
                    gsl::not_null<Cache*> cache,
                    hydro::Tags::Temperature<DataType> /*meta*/) const;
    void operator()(gsl::not_null<tnsr::II<DataType, 3>*> inv_spatial_metric,
                    gsl::not_null<Cache*> cache,
                    gr::Tags::InverseSpatialMetric<DataType, 3> /*meta*/) const;
    void operator()(
        gsl::not_null<tnsr::I<DataType, 3>*>
            lorentz_factor_times_spatial_velocity,
        gsl::not_null<Cache*> cache,
        hydro::Tags::LorentzFactorTimesSpatialVelocity<DataType, 3> /*meta*/)
        const;
    void operator()(gsl::not_null<tnsr::I<DataType, 3>*> spatial_velocity,
                    gsl::not_null<Cache*> cache,
                    hydro::Tags::SpatialVelocity<DataType, 3> /*meta*/) const;
    void operator()(gsl::not_null<Scalar<DataType>*> lorentz_factor,
                    gsl::not_null<Cache*> cache,
                    hydro::Tags::LorentzFactor<DataType> /*meta*/) const;
    void operator()(gsl::not_null<Scalar<DataType>*> electron_fraction,
                    gsl::not_null<Cache*> cache,
                    hydro::Tags::ElectronFraction<DataType> /*meta*/) const;
    void operator()(gsl::not_null<tnsr::I<DataType, 3>*> magnetic_field,
                    gsl::not_null<Cache*> cache,
                    hydro::Tags::MagneticField<DataType, 3> /*meta*/) const;
    void operator()(
        gsl::not_null<Scalar<DataType>*> div_cleaning_field,
        gsl::not_null<Cache*> cache,
        hydro::Tags::DivergenceCleaningField<DataType> /*meta*/) const;
  };

  std::string data_directory_{};
  std::unique_ptr<equation_of_state_type> equation_of_state_{nullptr};
  double density_cutoff_ = std::numeric_limits<double>::signaling_NaN();
  double atmosphere_density_ = std::numeric_limits<double>::signaling_NaN();
  std::optional<double> electron_fraction_ = std::nullopt;

  std::unique_ptr<spec::Exporter> spec_exporter_{nullptr};
};

}  // namespace grmhd::AnalyticData
