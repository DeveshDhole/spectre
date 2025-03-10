// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <algorithm>
#include <array>
#include <limits>
#include <memory>
#include <tuple>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Framework/CheckWithRandomValues.hpp"
#include "Framework/SetupLocalPythonEnvironment.hpp"
#include "Framework/TestCreation.hpp"
#include "Framework/TestHelpers.hpp"
#include "PointwiseFunctions/AnalyticData/GrMhd/MagneticRotor.hpp"
#include "PointwiseFunctions/Hydro/Tags.hpp"
#include "PointwiseFunctions/InitialDataUtilities/InitialData.hpp"
#include "PointwiseFunctions/InitialDataUtilities/Tags/InitialData.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace {

struct MagneticRotorProxy : grmhd::AnalyticData::MagneticRotor {
  using grmhd::AnalyticData::MagneticRotor::MagneticRotor;

  template <typename DataType>
  using hydro_variables_tags =
      tmpl::list<hydro::Tags::RestMassDensity<DataType>,
                 hydro::Tags::SpatialVelocity<DataType, 3>,
                 hydro::Tags::SpecificInternalEnergy<DataType>,
                 hydro::Tags::Pressure<DataType>,
                 hydro::Tags::LorentzFactor<DataType>,
                 hydro::Tags::SpecificEnthalpy<DataType>>;

  template <typename DataType>
  using grmhd_variables_tags =
      tmpl::push_back<hydro_variables_tags<DataType>,
                      hydro::Tags::MagneticField<DataType, 3>,
                      hydro::Tags::DivergenceCleaningField<DataType>>;

  template <typename DataType>
  tuples::tagged_tuple_from_typelist<hydro_variables_tags<DataType>>
  hydro_variables(const tnsr::I<DataType, 3, Frame::Inertial>& x) const {
    return variables(x, hydro_variables_tags<DataType>{});
  }

  template <typename DataType>
  tuples::tagged_tuple_from_typelist<grmhd_variables_tags<DataType>>
  grmhd_variables(const tnsr::I<DataType, 3, Frame::Inertial>& x) const {
    return variables(x, grmhd_variables_tags<DataType>{});
  }
};

void test_create_from_options() {
  register_classes_with_charm<grmhd::AnalyticData::MagneticRotor>();
  const std::unique_ptr<evolution::initial_data::InitialData> option_solution =
      TestHelpers::test_option_tag_factory_creation<
          evolution::initial_data::OptionTags::InitialData,
          grmhd::AnalyticData::MagneticRotor>(
          "MagneticRotor:\n"
          "  RotorRadius: 0.1\n"
          "  RotorDensity: 10.0\n"
          "  BackgroundDensity: 1.0\n"
          "  Pressure: 1.0\n"
          "  AngularVelocity: 9.95\n"
          "  MagneticField: [3.5449077018, 0.0, 0.0]\n"
          "  AdiabaticIndex: 1.6666666666666666")
          ->get_clone();
  const auto deserialized_option_solution =
      serialize_and_deserialize(option_solution);
  const auto& magnetic_rotor =
      dynamic_cast<const grmhd::AnalyticData::MagneticRotor&>(
          *deserialized_option_solution);

  CHECK(magnetic_rotor == grmhd::AnalyticData::MagneticRotor(
                              0.1, 10.0, 1.0, 1.0, 9.95,
                              std::array<double, 3>{{3.5449077018, 0.0, 0.0}},
                              1.6666666666666666));
}

void test_move() {
  grmhd::AnalyticData::MagneticRotor magnetic_rotor(
      0.1, 10.0, 1.0, 1.0, 9.95,
      std::array<double, 3>{{3.5449077018, 0.0, 0.0}}, 1.6666666666666666);
  grmhd::AnalyticData::MagneticRotor magnetic_rotor_copy(
      0.1, 10.0, 1.0, 1.0, 9.95,
      std::array<double, 3>{{3.5449077018, 0.0, 0.0}}, 1.6666666666666666);
  test_move_semantics(std::move(magnetic_rotor),
                      magnetic_rotor_copy);  //  NOLINT
}

void test_serialize() {
  grmhd::AnalyticData::MagneticRotor magnetic_rotor(
      0.1, 10.0, 1.0, 1.0, 9.95,
      std::array<double, 3>{{3.5449077018, 0.0, 0.0}}, 1.6666666666666666);
  test_serialization(magnetic_rotor);
}

template <typename DataType>
void test_variables(const DataType& used_for_size) {
  const double rotor_radius = 0.1;
  const double rotor_density = 10.0;
  const double background_density = 1.0;
  const double pressure = 1.0;
  const double angular_velocity = 9.95;
  const std::array<double, 3> magnetic_field{{3.5449077018, 0.0, 0.0}};
  const double adiabatic_index = 1.6666666666666666;

  const auto member_variables =
      std::make_tuple(rotor_radius, rotor_density, background_density, pressure,
                      angular_velocity, magnetic_field, adiabatic_index);

  MagneticRotorProxy magnetic_rotor(
      rotor_radius, rotor_density, background_density, pressure,
      angular_velocity, magnetic_field, adiabatic_index);

  pypp::check_with_random_values<1>(
      &MagneticRotorProxy::hydro_variables<DataType>, magnetic_rotor,
      "MagneticRotor",
      {"rest_mass_density", "spatial_velocity", "specific_internal_energy",
       "compute_pressure", "lorentz_factor", "specific_enthalpy"},
      {{{-1.0, 1.0}}}, member_variables, used_for_size);

  pypp::check_with_random_values<1>(
      &MagneticRotorProxy::grmhd_variables<DataType>, magnetic_rotor,
      "MagneticRotor",
      {"rest_mass_density", "spatial_velocity", "specific_internal_energy",
       "compute_pressure", "lorentz_factor", "specific_enthalpy",
       "magnetic_field", "divergence_cleaning_field"},
      {{{-1.0, 1.0}}}, member_variables, used_for_size);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.PointwiseFunctions.AnalyticData.GrMhd.MagneticRotor",
                  "[Unit][PointwiseFunctions]") {
  pypp::SetupLocalPythonEnvironment local_python_env{
      "PointwiseFunctions/AnalyticData/GrMhd"};

  test_create_from_options();
  test_serialize();
  test_move();

  test_variables(std::numeric_limits<double>::signaling_NaN());
  test_variables(DataVector(5));

  CHECK_THROWS_WITH(
      (grmhd::AnalyticData::MagneticRotor(
          0.2, 10.0, 1.0, 1.0, -9.95,
          std::array<double, 3>{{3.5449077018, 0.0, 0.0}}, 1.6666666666666666)),
      Catch::Matchers::ContainsSubstring(
          "MagneticRotor expects RotorRadius * | AngularVelocity | < 1"));
}
