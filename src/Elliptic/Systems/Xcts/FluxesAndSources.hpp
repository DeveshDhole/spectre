// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/Tensor/Tensor.hpp"
#include "Elliptic/Systems/Xcts/Geometry.hpp"
#include "Elliptic/Systems/Xcts/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/Divergence.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags/Conformal.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
class DataVector;
/// \endcond

namespace Xcts {

/// Indicates a subset of the XCTS equations
enum class Equations {
  /// Only the Hamiltonian constraint, solved for \f$\psi\f$
  Hamiltonian,
  /// Both the Hamiltonian constraint and the lapse equation, solved for
  /// \f$\psi\f$ and \f$\alpha\psi\f$
  HamiltonianAndLapse,
  /// The full XCTS equations, solved for \f$\psi\f$, \f$\alpha\psi\f$ and
  /// \f$\beta_\mathrm{excess}\f$
  HamiltonianLapseAndShift
};

/// The fluxes \f$F^i\f$ for the first-order formulation of the XCTS equations.
///
/// \see Xcts::FirstOrderSystem
template <Equations EnabledEquations, Geometry ConformalGeometry>
struct Fluxes;

/// \cond
template <>
struct Fluxes<Equations::Hamiltonian, Geometry::FlatCartesian> {
  using argument_tags = tmpl::list<>;
  using volume_tags = tmpl::list<>;
  using const_global_cache_tags = tmpl::list<>;
  static constexpr bool is_trivial = true;
  static constexpr bool is_discontinuous = false;
  static void apply(
      gsl::not_null<tnsr::I<DataVector, 3>*> flux_for_conformal_factor,
      const Scalar<DataVector>& conformal_factor_minus_one,
      const tnsr::i<DataVector, 3>& conformal_factor_gradient);
  static void apply(
      gsl::not_null<tnsr::I<DataVector, 3>*> flux_for_conformal_factor,
      const tnsr::i<DataVector, 3>& face_normal,
      const tnsr::I<DataVector, 3>& face_normal_vector,
      const Scalar<DataVector>& conformal_factor_minus_one);
};

template <>
struct Fluxes<Equations::Hamiltonian, Geometry::Curved> {
  using argument_tags =
      tmpl::list<Tags::InverseConformalMetric<DataVector, 3, Frame::Inertial>>;
  using volume_tags = tmpl::list<>;
  using const_global_cache_tags = tmpl::list<>;
  static constexpr bool is_trivial = true;
  static constexpr bool is_discontinuous = false;
  static void apply(
      gsl::not_null<tnsr::I<DataVector, 3>*> flux_for_conformal_factor,
      const tnsr::II<DataVector, 3>& inv_conformal_metric,
      const Scalar<DataVector>& conformal_factor_minus_one,
      const tnsr::i<DataVector, 3>& conformal_factor_gradient);
  static void apply(
      gsl::not_null<tnsr::I<DataVector, 3>*> flux_for_conformal_factor,
      const tnsr::II<DataVector, 3>& inv_conformal_metric,
      const tnsr::i<DataVector, 3>& face_normal,
      const tnsr::I<DataVector, 3>& face_normal_vector,
      const Scalar<DataVector>& conformal_factor_minus_one);
};

template <>
struct Fluxes<Equations::HamiltonianAndLapse, Geometry::FlatCartesian> {
  using argument_tags = tmpl::list<>;
  using volume_tags = tmpl::list<>;
  using const_global_cache_tags = tmpl::list<>;
  static constexpr bool is_trivial = true;
  static constexpr bool is_discontinuous = false;
  static void apply(
      gsl::not_null<tnsr::I<DataVector, 3>*> flux_for_conformal_factor,
      gsl::not_null<tnsr::I<DataVector, 3>*>
          flux_for_lapse_times_conformal_factor,
      const Scalar<DataVector>& conformal_factor_minus_one,
      const Scalar<DataVector>& lapse_times_conformal_factor_minus_one,
      const tnsr::i<DataVector, 3>& conformal_factor_gradient,
      const tnsr::i<DataVector, 3>& lapse_times_conformal_factor_gradient);
  static void apply(
      gsl::not_null<tnsr::I<DataVector, 3>*> flux_for_conformal_factor,
      gsl::not_null<tnsr::I<DataVector, 3>*>
          flux_for_lapse_times_conformal_factor,
      const tnsr::i<DataVector, 3>& face_normal,
      const tnsr::I<DataVector, 3>& face_normal_vector,
      const Scalar<DataVector>& conformal_factor_minus_one,
      const Scalar<DataVector>& lapse_times_conformal_factor_minus_one);
};

template <>
struct Fluxes<Equations::HamiltonianAndLapse, Geometry::Curved> {
  using argument_tags =
      tmpl::list<Tags::InverseConformalMetric<DataVector, 3, Frame::Inertial>>;
  using volume_tags = tmpl::list<>;
  using const_global_cache_tags = tmpl::list<>;
  static constexpr bool is_trivial = true;
  static constexpr bool is_discontinuous = false;
  static void apply(
      gsl::not_null<tnsr::I<DataVector, 3>*> flux_for_conformal_factor,
      gsl::not_null<tnsr::I<DataVector, 3>*>
          flux_for_lapse_times_conformal_factor,
      const tnsr::II<DataVector, 3>& inv_conformal_metric,
      const Scalar<DataVector>& conformal_factor_minus_one,
      const Scalar<DataVector>& lapse_times_conformal_factor_minus_one,
      const tnsr::i<DataVector, 3>& conformal_factor_gradient,
      const tnsr::i<DataVector, 3>& lapse_times_conformal_factor_gradient);
  static void apply(
      gsl::not_null<tnsr::I<DataVector, 3>*> flux_for_conformal_factor,
      gsl::not_null<tnsr::I<DataVector, 3>*>
          flux_for_lapse_times_conformal_factor,
      const tnsr::II<DataVector, 3>& inv_conformal_metric,
      const tnsr::i<DataVector, 3>& face_normal,
      const tnsr::I<DataVector, 3>& face_normal_vector,
      const Scalar<DataVector>& conformal_factor_minus_one,
      const Scalar<DataVector>& lapse_times_conformal_factor_minus_one);
};

template <>
struct Fluxes<Equations::HamiltonianLapseAndShift, Geometry::FlatCartesian> {
  using argument_tags = tmpl::list<>;
  using volume_tags = tmpl::list<>;
  using const_global_cache_tags = tmpl::list<>;
  static constexpr bool is_trivial = false;
  static constexpr bool is_discontinuous = false;
  static void apply(
      gsl::not_null<tnsr::I<DataVector, 3>*> flux_for_conformal_factor,
      gsl::not_null<tnsr::I<DataVector, 3>*>
          flux_for_lapse_times_conformal_factor,
      gsl::not_null<tnsr::II<DataVector, 3>*> longitudinal_shift_excess,
      const Scalar<DataVector>& conformal_factor_minus_one,
      const Scalar<DataVector>& lapse_times_conformal_factor_minus_one,
      const tnsr::I<DataVector, 3>& shift_excess,
      const tnsr::i<DataVector, 3>& conformal_factor_gradient,
      const tnsr::i<DataVector, 3>& lapse_times_conformal_factor_gradient,
      const tnsr::iJ<DataVector, 3>& deriv_shift_excess);
  static void apply(
      gsl::not_null<tnsr::I<DataVector, 3>*> flux_for_conformal_factor,
      gsl::not_null<tnsr::I<DataVector, 3>*>
          flux_for_lapse_times_conformal_factor,
      gsl::not_null<tnsr::II<DataVector, 3>*> longitudinal_shift_excess,
      const tnsr::i<DataVector, 3>& face_normal,
      const tnsr::I<DataVector, 3>& face_normal_vector,
      const Scalar<DataVector>& conformal_factor_minus_one,
      const Scalar<DataVector>& lapse_times_conformal_factor_minus_one,
      const tnsr::I<DataVector, 3>& shift_excess);
};

template <>
struct Fluxes<Equations::HamiltonianLapseAndShift, Geometry::Curved> {
  using argument_tags = tmpl::list<
      Tags::ConformalMetric<DataVector, 3, Frame::Inertial>,
      Tags::InverseConformalMetric<DataVector, 3, Frame::Inertial>,
      Tags::ConformalChristoffelSecondKind<DataVector, 3, Frame::Inertial>>;
  using volume_tags = tmpl::list<>;
  using const_global_cache_tags = tmpl::list<>;
  static constexpr bool is_trivial = false;
  static constexpr bool is_discontinuous = false;
  static void apply(
      gsl::not_null<tnsr::I<DataVector, 3>*> flux_for_conformal_factor,
      gsl::not_null<tnsr::I<DataVector, 3>*>
          flux_for_lapse_times_conformal_factor,
      gsl::not_null<tnsr::II<DataVector, 3>*> longitudinal_shift_excess,
      const tnsr::ii<DataVector, 3>& /*conformal_metric*/,
      const tnsr::II<DataVector, 3>& inv_conformal_metric,
      const tnsr::Ijj<DataVector, 3>& christoffel_second_kind,
      const Scalar<DataVector>& conformal_factor_minus_one,
      const Scalar<DataVector>& lapse_times_conformal_factor_minus_one,
      const tnsr::I<DataVector, 3>& shift_excess,
      const tnsr::i<DataVector, 3>& conformal_factor_gradient,
      const tnsr::i<DataVector, 3>& lapse_times_conformal_factor_gradient,
      const tnsr::iJ<DataVector, 3>& deriv_shift_excess);
  static void apply(
      gsl::not_null<tnsr::I<DataVector, 3>*> flux_for_conformal_factor,
      gsl::not_null<tnsr::I<DataVector, 3>*>
          flux_for_lapse_times_conformal_factor,
      gsl::not_null<tnsr::II<DataVector, 3>*> longitudinal_shift_excess,
      const tnsr::ii<DataVector, 3>& /*conformal_metric*/,
      const tnsr::II<DataVector, 3>& inv_conformal_metric,
      const tnsr::Ijj<DataVector, 3>& christoffel_second_kind,
      const tnsr::i<DataVector, 3>& face_normal,
      const tnsr::I<DataVector, 3>& face_normal_vector,
      const Scalar<DataVector>& conformal_factor_minus_one,
      const Scalar<DataVector>& lapse_times_conformal_factor_minus_one,
      const tnsr::I<DataVector, 3>& shift_excess);
};
/// \endcond

/// The sources \f$S\f$ for the first-order formulation of the XCTS equations.
///
/// \see Xcts::FirstOrderSystem
template <Equations EnabledEquations, Geometry ConformalGeometry,
          int ConformalMatterScale>
struct Sources;

/// \cond
template <int ConformalMatterScale>
struct Sources<Equations::Hamiltonian, Geometry::FlatCartesian,
               ConformalMatterScale> {
  using argument_tags = tmpl::list<
      gr::Tags::Conformal<gr::Tags::EnergyDensity<DataVector>,
                          ConformalMatterScale>,
      gr::Tags::TraceExtrinsicCurvature<DataVector>,
      Tags::LongitudinalShiftMinusDtConformalMetricOverLapseSquare<DataVector>>;
  static void apply(
      gsl::not_null<Scalar<DataVector>*> hamiltonian_constraint,
      const Scalar<DataVector>& conformal_energy_density,
      const Scalar<DataVector>& extrinsic_curvature_trace,
      const Scalar<DataVector>&
          longitudinal_shift_minus_dt_conformal_metric_over_lapse_square,
      const Scalar<DataVector>& conformal_factor_minus_one,
      const tnsr::I<DataVector, 3>& conformal_factor_flux);
};

template <int ConformalMatterScale>
struct Sources<Equations::Hamiltonian, Geometry::Curved, ConformalMatterScale> {
  using argument_tags = tmpl::push_back<
      typename Sources<Equations::Hamiltonian, Geometry::FlatCartesian,
                       ConformalMatterScale>::argument_tags,
      Tags::ConformalChristoffelContracted<DataVector, 3, Frame::Inertial>,
      Tags::ConformalRicciScalar<DataVector>>;
  using const_global_cache_tags = tmpl::list<>;
  static void apply(
      gsl::not_null<Scalar<DataVector>*> hamiltonian_constraint,
      const Scalar<DataVector>& conformal_energy_density,
      const Scalar<DataVector>& extrinsic_curvature_trace,
      const Scalar<DataVector>&
          longitudinal_shift_minus_dt_conformal_metric_over_lapse_square,
      const tnsr::i<DataVector, 3>& conformal_christoffel_contracted,
      const Scalar<DataVector>& conformal_ricci_scalar,
      const Scalar<DataVector>& conformal_factor_minus_one,
      const tnsr::I<DataVector, 3>& conformal_factor_flux);
};

template <int ConformalMatterScale>
struct Sources<Equations::HamiltonianAndLapse, Geometry::FlatCartesian,
               ConformalMatterScale> {
  using argument_tags = tmpl::list<
      gr::Tags::Conformal<gr::Tags::EnergyDensity<DataVector>,
                          ConformalMatterScale>,
      gr::Tags::Conformal<gr::Tags::StressTrace<DataVector>,
                          ConformalMatterScale>,
      gr::Tags::TraceExtrinsicCurvature<DataVector>,
      ::Tags::dt<gr::Tags::TraceExtrinsicCurvature<DataVector>>,
      Tags::LongitudinalShiftMinusDtConformalMetricSquare<DataVector>,
      Tags::ShiftDotDerivExtrinsicCurvatureTrace<DataVector>>;
  using const_global_cache_tags = tmpl::list<>;
  static void apply(
      gsl::not_null<Scalar<DataVector>*> hamiltonian_constraint,
      gsl::not_null<Scalar<DataVector>*> lapse_equation,
      const Scalar<DataVector>& conformal_energy_density,
      const Scalar<DataVector>& conformal_stress_trace,
      const Scalar<DataVector>& extrinsic_curvature_trace,
      const Scalar<DataVector>& dt_extrinsic_curvature_trace,
      const Scalar<DataVector>&
          longitudinal_shift_minus_dt_conformal_metric_square,
      const Scalar<DataVector>& shift_dot_deriv_extrinsic_curvature_trace,
      const Scalar<DataVector>& conformal_factor_minus_one,
      const Scalar<DataVector>& lapse_times_conformal_factor_minus_one,
      const tnsr::I<DataVector, 3>& conformal_factor_flux,
      const tnsr::I<DataVector, 3>& lapse_times_conformal_factor_flux);
};

template <int ConformalMatterScale>
struct Sources<Equations::HamiltonianAndLapse, Geometry::Curved,
               ConformalMatterScale> {
  using argument_tags = tmpl::push_back<
      typename Sources<Equations::HamiltonianAndLapse, Geometry::FlatCartesian,
                       ConformalMatterScale>::argument_tags,
      Tags::ConformalChristoffelContracted<DataVector, 3, Frame::Inertial>,
      Tags::ConformalRicciScalar<DataVector>>;
  using const_global_cache_tags = tmpl::list<>;
  static void apply(
      gsl::not_null<Scalar<DataVector>*> hamiltonian_constraint,
      gsl::not_null<Scalar<DataVector>*> lapse_equation,
      const Scalar<DataVector>& conformal_energy_density,
      const Scalar<DataVector>& conformal_stress_trace,
      const Scalar<DataVector>& extrinsic_curvature_trace,
      const Scalar<DataVector>& dt_extrinsic_curvature_trace,
      const Scalar<DataVector>&
          longitudinal_shift_minus_dt_conformal_metric_square,
      const Scalar<DataVector>& shift_dot_deriv_extrinsic_curvature_trace,
      const tnsr::i<DataVector, 3>& conformal_christoffel_contracted,
      const Scalar<DataVector>& conformal_ricci_scalar,
      const Scalar<DataVector>& conformal_factor_minus_one,
      const Scalar<DataVector>& lapse_times_conformal_factor_minus_one,
      const tnsr::I<DataVector, 3>& conformal_factor_flux,
      const tnsr::I<DataVector, 3>& lapse_times_conformal_factor_flux);
};

template <int ConformalMatterScale>
struct Sources<Equations::HamiltonianLapseAndShift, Geometry::FlatCartesian,
               ConformalMatterScale> {
  using argument_tags = tmpl::list<
      gr::Tags::Conformal<gr::Tags::EnergyDensity<DataVector>,
                          ConformalMatterScale>,
      gr::Tags::Conformal<gr::Tags::StressTrace<DataVector>,
                          ConformalMatterScale>,
      gr::Tags::Conformal<gr::Tags::MomentumDensity<DataVector, 3>,
                          ConformalMatterScale>,
      gr::Tags::TraceExtrinsicCurvature<DataVector>,
      ::Tags::dt<gr::Tags::TraceExtrinsicCurvature<DataVector>>,
      ::Tags::deriv<gr::Tags::TraceExtrinsicCurvature<DataVector>,
                    tmpl::size_t<3>, Frame::Inertial>,
      Tags::ShiftBackground<DataVector, 3, Frame::Inertial>,
      Tags::LongitudinalShiftBackgroundMinusDtConformalMetric<DataVector, 3,
                                                              Frame::Inertial>,
      ::Tags::div<Tags::LongitudinalShiftBackgroundMinusDtConformalMetric<
          DataVector, 3, Frame::Inertial>>>;
  using const_global_cache_tags = tmpl::list<>;
  static void apply(
      gsl::not_null<Scalar<DataVector>*> hamiltonian_constraint,
      gsl::not_null<Scalar<DataVector>*> lapse_equation,
      gsl::not_null<tnsr::I<DataVector, 3>*> momentum_constraint,
      const Scalar<DataVector>& conformal_energy_density,
      const Scalar<DataVector>& conformal_stress_trace,
      const tnsr::I<DataVector, 3>& conformal_momentum_density,
      const Scalar<DataVector>& extrinsic_curvature_trace,
      const Scalar<DataVector>& dt_extrinsic_curvature_trace,
      const tnsr::i<DataVector, 3>& extrinsic_curvature_trace_gradient,
      const tnsr::I<DataVector, 3>& shift_background,
      const tnsr::II<DataVector, 3>&
          longitudinal_shift_background_minus_dt_conformal_metric,
      const tnsr::I<DataVector, 3>&
          div_longitudinal_shift_background_minus_dt_conformal_metric,
      const Scalar<DataVector>& conformal_factor_minus_one,
      const Scalar<DataVector>& lapse_times_conformal_factor_minus_one,
      const tnsr::I<DataVector, 3>& shift_excess,
      const tnsr::I<DataVector, 3>& conformal_factor_flux,
      const tnsr::I<DataVector, 3>& lapse_times_conformal_factor_flux,
      const tnsr::II<DataVector, 3>& longitudinal_shift_excess);
};

template <int ConformalMatterScale>
struct Sources<Equations::HamiltonianLapseAndShift, Geometry::Curved,
               ConformalMatterScale> {
  using argument_tags = tmpl::push_back<
      typename Sources<Equations::HamiltonianLapseAndShift,
                       Geometry::FlatCartesian,
                       ConformalMatterScale>::argument_tags,
      Tags::ConformalMetric<DataVector, 3, Frame::Inertial>,
      Tags::InverseConformalMetric<DataVector, 3, Frame::Inertial>,
      Tags::ConformalChristoffelFirstKind<DataVector, 3, Frame::Inertial>,
      Tags::ConformalChristoffelSecondKind<DataVector, 3, Frame::Inertial>,
      Tags::ConformalChristoffelContracted<DataVector, 3, Frame::Inertial>,
      Tags::ConformalRicciScalar<DataVector>>;
  using const_global_cache_tags = tmpl::list<>;
  static void apply(
      gsl::not_null<Scalar<DataVector>*> hamiltonian_constraint,
      gsl::not_null<Scalar<DataVector>*> lapse_equation,
      gsl::not_null<tnsr::I<DataVector, 3>*> momentum_constraint,
      const Scalar<DataVector>& conformal_energy_density,
      const Scalar<DataVector>& conformal_stress_trace,
      const tnsr::I<DataVector, 3>& conformal_momentum_density,
      const Scalar<DataVector>& extrinsic_curvature_trace,
      const Scalar<DataVector>& dt_extrinsic_curvature_trace,
      const tnsr::i<DataVector, 3>& extrinsic_curvature_trace_gradient,
      const tnsr::I<DataVector, 3>& shift_background,
      const tnsr::II<DataVector, 3>&
          longitudinal_shift_background_minus_dt_conformal_metric,
      const tnsr::I<DataVector, 3>&
          div_longitudinal_shift_background_minus_dt_conformal_metric,
      const tnsr::ii<DataVector, 3>& conformal_metric,
      const tnsr::II<DataVector, 3>& inv_conformal_metric,
      const tnsr::ijj<DataVector, 3>& /*conformal_christoffel_first_kind*/,
      const tnsr::Ijj<DataVector, 3>& conformal_christoffel_second_kind,
      const tnsr::i<DataVector, 3>& conformal_christoffel_contracted,
      const Scalar<DataVector>& conformal_ricci_scalar,
      const Scalar<DataVector>& conformal_factor_minus_one,
      const Scalar<DataVector>& lapse_times_conformal_factor_minus_one,
      const tnsr::I<DataVector, 3>& shift_excess,
      const tnsr::I<DataVector, 3>& conformal_factor_flux,
      const tnsr::I<DataVector, 3>& lapse_times_conformal_factor_flux,
      const tnsr::II<DataVector, 3>& longitudinal_shift_excess);
};
/// \endcond

/// The linearization of the sources \f$S\f$ for the first-order formulation of
/// the XCTS equations.
///
/// \see Xcts::FirstOrderSystem
template <Equations EnabledEquations, Geometry ConformalGeometry,
          int ConformalMatterScale>
struct LinearizedSources;

/// \cond
template <int ConformalMatterScale>
struct LinearizedSources<Equations::Hamiltonian, Geometry::FlatCartesian,
                         ConformalMatterScale> {
  using argument_tags = tmpl::push_back<
      typename Sources<Equations::Hamiltonian, Geometry::FlatCartesian,
                       ConformalMatterScale>::argument_tags,
      Tags::ConformalFactorMinusOne<DataVector>>;
  using const_global_cache_tags = tmpl::list<>;
  static void apply(
      gsl::not_null<Scalar<DataVector>*> linearized_hamiltonian_constraint,
      const Scalar<DataVector>& conformal_energy_density,
      const Scalar<DataVector>& extrinsic_curvature_trace,
      const Scalar<DataVector>&
          longitudinal_shift_minus_dt_conformal_metric_over_lapse_square,
      const Scalar<DataVector>& conformal_factor_minus_one,
      const Scalar<DataVector>& conformal_factor_correction,
      const tnsr::I<DataVector, 3>&
      /*conformal_factor_flux_correction*/);
};

template <int ConformalMatterScale>
struct LinearizedSources<Equations::Hamiltonian, Geometry::Curved,
                         ConformalMatterScale> {
  using argument_tags =
      tmpl::push_back<typename Sources<Equations::Hamiltonian, Geometry::Curved,
                                       ConformalMatterScale>::argument_tags,
                      Tags::ConformalFactorMinusOne<DataVector>>;
  using const_global_cache_tags = tmpl::list<>;
  static void apply(
      gsl::not_null<Scalar<DataVector>*> linearized_hamiltonian_constraint,
      const Scalar<DataVector>& conformal_energy_density,
      const Scalar<DataVector>& extrinsic_curvature_trace,
      const Scalar<DataVector>&
          longitudinal_shift_minus_dt_conformal_metric_over_lapse_square,
      const tnsr::i<DataVector, 3>& conformal_christoffel_contracted,
      const Scalar<DataVector>& conformal_ricci_scalar,
      const Scalar<DataVector>& conformal_factor_minus_one,
      const Scalar<DataVector>& conformal_factor_correction,
      const tnsr::I<DataVector, 3>& conformal_factor_flux_correction);
};

template <int ConformalMatterScale>
struct LinearizedSources<Equations::HamiltonianAndLapse,
                         Geometry::FlatCartesian, ConformalMatterScale> {
  using argument_tags = tmpl::push_back<
      typename Sources<Equations::HamiltonianAndLapse, Geometry::FlatCartesian,
                       ConformalMatterScale>::argument_tags,
      Tags::ConformalFactorMinusOne<DataVector>,
      Tags::LapseTimesConformalFactorMinusOne<DataVector>>;
  using const_global_cache_tags = tmpl::list<>;
  static void apply(
      gsl::not_null<Scalar<DataVector>*> linearized_hamiltonian_constraint,
      gsl::not_null<Scalar<DataVector>*> linearized_lapse_equation,
      const Scalar<DataVector>& conformal_energy_density,
      const Scalar<DataVector>& conformal_stress_trace,
      const Scalar<DataVector>& extrinsic_curvature_trace,
      const Scalar<DataVector>& dt_extrinsic_curvature_trace,
      const Scalar<DataVector>&
          longitudinal_shift_minus_dt_conformal_metric_square,
      const Scalar<DataVector>& shift_dot_deriv_extrinsic_curvature_trace,
      const Scalar<DataVector>& conformal_factor_minus_one,
      const Scalar<DataVector>& lapse_times_conformal_factor_minus_one,
      const Scalar<DataVector>& conformal_factor_correction,
      const Scalar<DataVector>& lapse_times_conformal_factor_correction,
      const tnsr::I<DataVector, 3>& /*conformal_factor_flux_correction*/,
      const tnsr::I<DataVector, 3>&
          lapse_times_conformal_factor_flux_correction);
};

template <int ConformalMatterScale>
struct LinearizedSources<Equations::HamiltonianAndLapse, Geometry::Curved,
                         ConformalMatterScale> {
  using argument_tags = tmpl::push_back<
      typename Sources<Equations::HamiltonianAndLapse, Geometry::Curved,
                       ConformalMatterScale>::argument_tags,
      Tags::ConformalFactorMinusOne<DataVector>,
      Tags::LapseTimesConformalFactorMinusOne<DataVector>>;
  using const_global_cache_tags = tmpl::list<>;
  static void apply(
      gsl::not_null<Scalar<DataVector>*> linearized_hamiltonian_constraint,
      gsl::not_null<Scalar<DataVector>*> linearized_lapse_equation,
      const Scalar<DataVector>& conformal_energy_density,
      const Scalar<DataVector>& conformal_stress_trace,
      const Scalar<DataVector>& extrinsic_curvature_trace,
      const Scalar<DataVector>& dt_extrinsic_curvature_trace,
      const Scalar<DataVector>&
          longitudinal_shift_minus_dt_conformal_metric_square,
      const Scalar<DataVector>& shift_dot_deriv_extrinsic_curvature_trace,
      const tnsr::i<DataVector, 3>& conformal_christoffel_contracted,
      const Scalar<DataVector>& conformal_ricci_scalar,
      const Scalar<DataVector>& conformal_factor_minus_one,
      const Scalar<DataVector>& lapse_times_conformal_factor_minus_one,
      const Scalar<DataVector>& conformal_factor_correction,
      const Scalar<DataVector>& lapse_times_conformal_factor_correction,
      const tnsr::I<DataVector, 3>& conformal_factor_flux_correction,
      const tnsr::I<DataVector, 3>&
          lapse_times_conformal_factor_flux_correction);
};

template <int ConformalMatterScale>
struct LinearizedSources<Equations::HamiltonianLapseAndShift,
                         Geometry::FlatCartesian, ConformalMatterScale> {
  using argument_tags = tmpl::push_back<
      typename Sources<Equations::HamiltonianLapseAndShift,
                       Geometry::FlatCartesian,
                       ConformalMatterScale>::argument_tags,
      Tags::ConformalFactorMinusOne<DataVector>,
      Tags::LapseTimesConformalFactorMinusOne<DataVector>,
      Tags::ShiftExcess<DataVector, 3, Frame::Inertial>,
      ::Tags::Flux<Tags::ConformalFactorMinusOne<DataVector>, tmpl::size_t<3>,
                   Frame::Inertial>,
      ::Tags::Flux<Tags::LapseTimesConformalFactorMinusOne<DataVector>,
                   tmpl::size_t<3>, Frame::Inertial>,
      Tags::LongitudinalShiftExcess<DataVector, 3, Frame::Inertial>>;
  using const_global_cache_tags = tmpl::list<>;
  static void apply(
      gsl::not_null<Scalar<DataVector>*> linearized_hamiltonian_constraint,
      gsl::not_null<Scalar<DataVector>*> linearized_lapse_equation,
      gsl::not_null<tnsr::I<DataVector, 3>*> linearized_momentum_constraint,
      const Scalar<DataVector>& conformal_energy_density,
      const Scalar<DataVector>& conformal_stress_trace,
      const tnsr::I<DataVector, 3>& conformal_momentum_density,
      const Scalar<DataVector>& extrinsic_curvature_trace,
      const Scalar<DataVector>& dt_extrinsic_curvature_trace,
      const tnsr::i<DataVector, 3>& extrinsic_curvature_trace_gradient,
      const tnsr::I<DataVector, 3>& shift_background,
      const tnsr::II<DataVector, 3>&
          longitudinal_shift_background_minus_dt_conformal_metric,
      const tnsr::I<DataVector, 3>&
          div_longitudinal_shift_background_minus_dt_conformal_metric,
      const Scalar<DataVector>& conformal_factor_minus_one,
      const Scalar<DataVector>& lapse_times_conformal_factor_minus_one,
      const tnsr::I<DataVector, 3>& shift_excess,
      const tnsr::I<DataVector, 3>& conformal_factor_flux,
      const tnsr::I<DataVector, 3>& lapse_times_conformal_factor_flux,
      const tnsr::II<DataVector, 3>& longitudinal_shift_excess,
      const Scalar<DataVector>& conformal_factor_correction,
      const Scalar<DataVector>& lapse_times_conformal_factor_correction,
      const tnsr::I<DataVector, 3>& shift_excess_correction,
      const tnsr::I<DataVector, 3>& conformal_factor_flux_correction,
      const tnsr::I<DataVector, 3>&
          lapse_times_conformal_factor_flux_correction,
      const tnsr::II<DataVector, 3>& longitudinal_shift_excess_correction);
};

template <int ConformalMatterScale>
struct LinearizedSources<Equations::HamiltonianLapseAndShift, Geometry::Curved,
                         ConformalMatterScale> {
  using argument_tags = tmpl::push_back<
      typename Sources<Equations::HamiltonianLapseAndShift, Geometry::Curved,
                       ConformalMatterScale>::argument_tags,
      Tags::ConformalFactorMinusOne<DataVector>,
      Tags::LapseTimesConformalFactorMinusOne<DataVector>,
      Tags::ShiftExcess<DataVector, 3, Frame::Inertial>,
      ::Tags::Flux<Tags::ConformalFactorMinusOne<DataVector>, tmpl::size_t<3>,
                   Frame::Inertial>,
      ::Tags::Flux<Tags::LapseTimesConformalFactorMinusOne<DataVector>,
                   tmpl::size_t<3>, Frame::Inertial>,
      Tags::LongitudinalShiftExcess<DataVector, 3, Frame::Inertial>>;
  using const_global_cache_tags = tmpl::list<>;
  static void apply(
      gsl::not_null<Scalar<DataVector>*> linearized_hamiltonian_constraint,
      gsl::not_null<Scalar<DataVector>*> linearized_lapse_equation,
      gsl::not_null<tnsr::I<DataVector, 3>*> linearized_momentum_constraint,
      const Scalar<DataVector>& conformal_energy_density,
      const Scalar<DataVector>& conformal_stress_trace,
      const tnsr::I<DataVector, 3>& conformal_momentum_density,
      const Scalar<DataVector>& extrinsic_curvature_trace,
      const Scalar<DataVector>& dt_extrinsic_curvature_trace,
      const tnsr::i<DataVector, 3>& extrinsic_curvature_trace_gradient,
      const tnsr::I<DataVector, 3>& shift_background,
      const tnsr::II<DataVector, 3>&
          longitudinal_shift_background_minus_dt_conformal_metric,
      const tnsr::I<DataVector, 3>&
      /*div_longitudinal_shift_background_minus_dt_conformal_metric*/,
      const tnsr::ii<DataVector, 3>& conformal_metric,
      const tnsr::II<DataVector, 3>& inv_conformal_metric,
      const tnsr::ijj<DataVector, 3>& /*conformal_christoffel_first_kind*/,
      const tnsr::Ijj<DataVector, 3>& conformal_christoffel_second_kind,
      const tnsr::i<DataVector, 3>& conformal_christoffel_contracted,
      const Scalar<DataVector>& conformal_ricci_scalar,
      const Scalar<DataVector>& conformal_factor_minus_one,
      const Scalar<DataVector>& lapse_times_conformal_factor_minus_one,
      const tnsr::I<DataVector, 3>& shift_excess,
      const tnsr::I<DataVector, 3>& conformal_factor_flux,
      const tnsr::I<DataVector, 3>& lapse_times_conformal_factor_flux,
      const tnsr::II<DataVector, 3>& longitudinal_shift_excess,
      const Scalar<DataVector>& conformal_factor_correction,
      const Scalar<DataVector>& lapse_times_conformal_factor_correction,
      const tnsr::I<DataVector, 3>& shift_excess_correction,
      const tnsr::I<DataVector, 3>& conformal_factor_flux_correction,
      const tnsr::I<DataVector, 3>&
          lapse_times_conformal_factor_flux_correction,
      const tnsr::II<DataVector, 3>& longitudinal_shift_excess_correction);
};
/// \endcond

}  // namespace Xcts
