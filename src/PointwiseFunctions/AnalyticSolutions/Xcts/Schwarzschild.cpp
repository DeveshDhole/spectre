// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "PointwiseFunctions/AnalyticSolutions/Xcts/Schwarzschild.hpp"

#include <algorithm>
#include <cstddef>
#include <limits>
#include <ostream>
#include <pup.h>
#include <string>
#include <utility>

#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/EagerMath/Magnitude.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Elliptic/Systems/Xcts/Tags.hpp"
#include "NumericalAlgorithms/LinearOperators/PartialDerivatives.hpp"
#include "NumericalAlgorithms/RootFinding/TOMS748.hpp"
#include "Options/Options.hpp"
#include "Options/ParseOptions.hpp"
#include "PointwiseFunctions/AnalyticSolutions/Xcts/CommonVariables.tpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags/Conformal.hpp"
#include "PointwiseFunctions/Xcts/ExtrinsicCurvature.hpp"
#include "Utilities/ConstantExpressions.hpp"
#include "Utilities/ErrorHandling/Error.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeWithValue.hpp"
#include "Utilities/Simd/Simd.hpp"
#include "Utilities/TMPL.hpp"

namespace Xcts::Solutions {

std::ostream& operator<<(std::ostream& os,
                         const SchwarzschildCoordinates coords) {
  switch (coords) {
    case SchwarzschildCoordinates::Isotropic:
      return os << "Isotropic";
    case SchwarzschildCoordinates::PainleveGullstrand:
      return os << "PainleveGullstrand";
    case SchwarzschildCoordinates::KerrSchildIsotropic:
      return os << "KerrSchildIsotropic";
    case SchwarzschildCoordinates::MaximalIsotropic:
      return os << "MaximalIsotropic";
      // LCOV_EXCL_START
    default:
      ERROR("Unknown SchwarzschildCoordinates");
      // LCOV_EXCL_END
  }
}

}  // namespace Xcts::Solutions

template <>
Xcts::Solutions::SchwarzschildCoordinates
Options::create_from_yaml<Xcts::Solutions::SchwarzschildCoordinates>::create<
    void>(const Options::Option& options) {
  const auto type_read = options.parse_as<std::string>();
  if ("Isotropic" == type_read) {
    return Xcts::Solutions::SchwarzschildCoordinates::Isotropic;
  } else if ("PainleveGullstrand" == type_read) {
    return Xcts::Solutions::SchwarzschildCoordinates::PainleveGullstrand;
  } else if ("KerrSchildIsotropic" == type_read) {
    return Xcts::Solutions::SchwarzschildCoordinates::KerrSchildIsotropic;
  } else if ("MaximalIsotropic" == type_read) {
    return Xcts::Solutions::SchwarzschildCoordinates::MaximalIsotropic;
  }
  PARSE_ERROR(options.context(),
              "Failed to convert \""
                  << type_read
                  << "\" to Xcts::Solutions::SchwarzschildCoordinates. Must be "
                     "one of 'Isotropic', 'PainleveGullstrand', "
                     "'KerrSchildIsotropic', 'MaximalIsotropic'.");
}

namespace Xcts::Solutions {
namespace detail {

namespace {
template <typename DataType>
DataType kerr_schild_isotropic_radius_from_areal(const DataType& areal_radius,
                                                 const double mass) {
  // Eq. (7.34) in https://arxiv.org/abs/gr-qc/0510016
  const auto one_over_lapse = sqrt(1. + 2. * mass / areal_radius);
  return 0.25 * areal_radius * square(1. + one_over_lapse) *
         exp(2. - 2. * one_over_lapse);
}

template <typename DataType>
DataType kerr_schild_isotropic_radius_from_areal_deriv(
    const DataType& areal_radius, const double mass) {
  const auto isotropic_radius =
      kerr_schild_isotropic_radius_from_areal(areal_radius, mass);
  const auto one_over_lapse = sqrt(1. + 2. * mass / areal_radius);
  return isotropic_radius / areal_radius * one_over_lapse;
}

double kerr_schild_areal_radius_from_isotropic(const double isotropic_radius,
                                               const double mass) {
  return RootFinder::toms748(
      [&isotropic_radius, &mass](const double areal_radius) {
        return kerr_schild_isotropic_radius_from_areal(areal_radius, mass) -
               isotropic_radius;
      },
      isotropic_radius, isotropic_radius + mass, 1.0e-12, 1.0e-15);
}

DataVector kerr_schild_areal_radius_from_isotropic(
    const DataVector& isotropic_radius, const double mass) {
  return RootFinder::toms748<true>(
      [&isotropic_radius, &mass](const auto areal_radius, const size_t i) {
        if constexpr (simd::is_batch<
                          std::decay_t<decltype(areal_radius)>>::value) {
          return kerr_schild_isotropic_radius_from_areal(areal_radius, mass) -
                 simd::load_unaligned(&isotropic_radius[i]);
        } else {
          return kerr_schild_isotropic_radius_from_areal(areal_radius, mass) -
                 isotropic_radius[i];
        }
      },
      isotropic_radius, isotropic_radius + mass, 1.0e-12, 1.0e-15);
}

template <typename DataType>
DataType maximal_isotropic_radius_from_areal(const DataType& areal_radius,
                                             const double mass) {
  return .25 *
         (2. * areal_radius + mass +
          sqrt(4. * square(areal_radius) + 4. * areal_radius * mass +
               3. * square(mass))) *
         pow((4. + 3. * sqrt(2.)) * (2. * areal_radius - 3. * mass) /
                 (8. * areal_radius + 6. * mass +
                  3. * sqrt(8. * square(areal_radius) +
                            8. * areal_radius * mass + 6. * square(mass))),
             1. / sqrt(2.));
}

template <typename DataType>
DataType maximal_isotropic_radius_from_areal_deriv(
    const DataType& areal_radius, const DataType& isotropic_radius,
    const double mass) {
  // r_isotropic = (A/4) * B^(1/sqrt(2))
  const DataType S = sqrt(8. * square(areal_radius) + 8. * areal_radius * mass +
                          6. * square(mass));
  const double C = 4. + 3. * sqrt(2.);
  const DataType D = 8. * areal_radius + 6. * mass + 3. * S;
  const DataType E = 2. * areal_radius - 3. * mass;
  const DataType F = 2. * areal_radius + mass;
  const DataType A = F + S / sqrt(2.);
  const DataType B = C * E / D;
  const DataType dAdR = 2. + (4. / sqrt(2.)) * F / S;
  const DataType dBdR = C * (2. * D - E * (8. + 12. * F / S)) / square(D);

  return isotropic_radius * (dAdR / A + dBdR / (B * sqrt(2.)));
}

template <typename DataType>
DataType areal_radius_from_maximal_isotropic(const DataType& isotropic_radius,
                                             const double mass) {
  const auto residual = [&isotropic_radius, &mass](const auto areal_radius,
                                                   const size_t i = 0) {
    if constexpr (simd::is_batch<std::decay_t<decltype(areal_radius)>>::value) {
      return maximal_isotropic_radius_from_areal(areal_radius, mass) -
             simd::load_unaligned(&(get_element(isotropic_radius, i)));
    } else {
      return maximal_isotropic_radius_from_areal(areal_radius, mass) -
             get_element(isotropic_radius, i);
    }
  };
  const auto lower_bound =
      make_with_value<DataType>(isotropic_radius, 1.5 * mass);
  const auto upper_bound = make_with_value<DataType>(isotropic_radius, 1e9);
  return RootFinder::toms748<false>(residual, lower_bound, upper_bound, 1.0e-12,
                                    1.0e-15);
}

}  // namespace

SchwarzschildImpl::SchwarzschildImpl(
    const double mass, const SchwarzschildCoordinates coordinate_system)
    : mass_(mass), coordinate_system_(coordinate_system) {}

double SchwarzschildImpl::mass() const { return mass_; }

SchwarzschildCoordinates SchwarzschildImpl::coordinate_system() const {
  return coordinate_system_;
}

double SchwarzschildImpl::radius_at_horizon() const {
  switch (coordinate_system_) {
    case SchwarzschildCoordinates::Isotropic:
      return 0.5 * mass_;
    case SchwarzschildCoordinates::PainleveGullstrand:
      return 2. * mass_;
    case SchwarzschildCoordinates::KerrSchildIsotropic:
      return kerr_schild_isotropic_radius_from_areal(2. * mass_, mass_);
    case SchwarzschildCoordinates::MaximalIsotropic:
      return maximal_isotropic_radius_from_areal(2. * mass_, mass_);
      // LCOV_EXCL_START
    default:
      ERROR("Missing case for SchwarzschildCoordinates");
      // LCOV_EXCL_END
  }
}

void SchwarzschildImpl::pup(PUP::er& p) {
  p | mass_;
  p | coordinate_system_;
}

bool operator==(const SchwarzschildImpl& lhs, const SchwarzschildImpl& rhs) {
  return lhs.mass() == rhs.mass() and
         lhs.coordinate_system() == rhs.coordinate_system();
}

bool operator!=(const SchwarzschildImpl& lhs, const SchwarzschildImpl& rhs) {
  return not(lhs == rhs);
}

template <typename DataType>
void SchwarzschildVariables<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> radius,
    const gsl::not_null<Cache*> /*cache*/,
    detail::Tags::Radius<DataType> /*meta*/) const {
  magnitude(radius, x);
}

template <typename DataType>
void SchwarzschildVariables<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> areal_radius,
    const gsl::not_null<Cache*> cache,
    detail::Tags::ArealRadius<DataType> /*meta*/) const {
  const auto& isotropic_radius =
      cache->get_var(*this, detail::Tags::Radius<DataType>{});
  switch (coordinate_system) {
    case SchwarzschildCoordinates::KerrSchildIsotropic:
      get(*areal_radius) =
          kerr_schild_areal_radius_from_isotropic(get(isotropic_radius), mass);
      break;
    case SchwarzschildCoordinates::MaximalIsotropic:
      get(*areal_radius) =
          areal_radius_from_maximal_isotropic(get(isotropic_radius), mass);
      break;
    default:
      ERROR(
          "The areal radius is only needed for 'KerrSchildIsotropic' "
          "or 'MaximalIsotropic' coordinates.");
  }
}

template <typename DataType>
void SchwarzschildVariables<DataType>::operator()(
    const gsl::not_null<tnsr::ii<DataType, 3>*> conformal_metric,
    const gsl::not_null<Cache*> /*cache*/,
    Xcts::Tags::ConformalMetric<DataType, 3, Frame::Inertial> /*meta*/) const {
  get<0, 0>(*conformal_metric) = 1.;
  get<1, 1>(*conformal_metric) = 1.;
  get<2, 2>(*conformal_metric) = 1.;
  get<0, 1>(*conformal_metric) = 0.;
  get<0, 2>(*conformal_metric) = 0.;
  get<1, 2>(*conformal_metric) = 0.;
}

template <typename DataType>
void SchwarzschildVariables<DataType>::operator()(
    const gsl::not_null<tnsr::II<DataType, 3>*> inv_conformal_metric,
    const gsl::not_null<Cache*> /*cache*/,
    Xcts::Tags::InverseConformalMetric<DataType, 3, Frame::Inertial> /*meta*/)
    const {
  get<0, 0>(*inv_conformal_metric) = 1.;
  get<1, 1>(*inv_conformal_metric) = 1.;
  get<2, 2>(*inv_conformal_metric) = 1.;
  get<0, 1>(*inv_conformal_metric) = 0.;
  get<0, 2>(*inv_conformal_metric) = 0.;
  get<1, 2>(*inv_conformal_metric) = 0.;
}

template <typename DataType>
void SchwarzschildVariables<DataType>::operator()(
    const gsl::not_null<tnsr::ijj<DataType, 3>*> deriv_conformal_metric,
    const gsl::not_null<Cache*> /*cache*/,
    ::Tags::deriv<Xcts::Tags::ConformalMetric<DataType, 3, Frame::Inertial>,
                  tmpl::size_t<3>, Frame::Inertial> /*meta*/) const {
  std::fill(deriv_conformal_metric->begin(), deriv_conformal_metric->end(), 0.);
}

template <typename DataType>
void SchwarzschildVariables<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> trace_extrinsic_curvature,
    const gsl::not_null<Cache*> cache,
    gr::Tags::TraceExtrinsicCurvature<DataType> /*meta*/) const {
  switch (coordinate_system) {
    case SchwarzschildCoordinates::Isotropic: {
      get(*trace_extrinsic_curvature) = 0.;
      break;
    }
    case SchwarzschildCoordinates::PainleveGullstrand: {
      const auto& r =
          get(cache->get_var(*this, detail::Tags::Radius<DataType>{}));
      get(*trace_extrinsic_curvature) = 1.5 * sqrt(2. * mass) / pow(r, 1.5);
      break;
    }
    case SchwarzschildCoordinates::KerrSchildIsotropic: {
      const auto& r =
          get(cache->get_var(*this, detail::Tags::ArealRadius<DataType>{}));
      const auto& lapse =
          get(cache->get_var(*this, gr::Tags::Lapse<DataType>{}));
      get(*trace_extrinsic_curvature) =
          2. * mass * cube(lapse) / square(r) * (1. + 3. * mass / r);
      break;
    }
    case SchwarzschildCoordinates::MaximalIsotropic: {
      get(*trace_extrinsic_curvature) = 0.;
      break;
    }
      // LCOV_EXCL_START
    default:
      ERROR("Missing case for SchwarzschildCoordinates");
      // LCOV_EXCL_END
  }
}

template <typename DataType>
void SchwarzschildVariables<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> dt_trace_extrinsic_curvature,
    const gsl::not_null<Cache*> /*cache*/,
    ::Tags::dt<gr::Tags::TraceExtrinsicCurvature<DataType>> /*meta*/) const {
  get(*dt_trace_extrinsic_curvature) = 0.;
}

template <typename DataType>
void SchwarzschildVariables<DataType>::operator()(
    const gsl::not_null<tnsr::i<DataType, 3>*>
        trace_extrinsic_curvature_gradient,
    const gsl::not_null<Cache*> cache,
    ::Tags::deriv<gr::Tags::TraceExtrinsicCurvature<DataType>, tmpl::size_t<3>,
                  Frame::Inertial> /*meta*/) const {
  switch (coordinate_system) {
    case SchwarzschildCoordinates::Isotropic: {
      std::fill(trace_extrinsic_curvature_gradient->begin(),
                trace_extrinsic_curvature_gradient->end(), 0.);
      break;
    }
    case SchwarzschildCoordinates::PainleveGullstrand: {
      const auto& r =
          get(cache->get_var(*this, detail::Tags::Radius<DataType>{}));
      const DataType isotropic_prefactor =
          -2.25 * sqrt(2. * mass) / pow(r, 3.5);
      get<0>(*trace_extrinsic_curvature_gradient) =
          isotropic_prefactor * get<0>(x);
      get<1>(*trace_extrinsic_curvature_gradient) =
          isotropic_prefactor * get<1>(x);
      get<2>(*trace_extrinsic_curvature_gradient) =
          isotropic_prefactor * get<2>(x);
      break;
    }
    case SchwarzschildCoordinates::KerrSchildIsotropic: {
      const auto& rbar =
          get(cache->get_var(*this, detail::Tags::Radius<DataType>{}));
      const auto& r =
          get(cache->get_var(*this, detail::Tags::ArealRadius<DataType>{}));
      const auto& lapse =
          get(cache->get_var(*this, gr::Tags::Lapse<DataType>{}));
      const auto& K = get(
          cache->get_var(*this, gr::Tags::TraceExtrinsicCurvature<DataType>{}));
      const DataType isotropic_prefactor =
          K * lapse *
          (3. * mass * square(lapse) / r - 2. - 3. * mass / (r + 3. * mass)) /
          square(rbar);
      get<0>(*trace_extrinsic_curvature_gradient) =
          isotropic_prefactor * get<0>(x);
      get<1>(*trace_extrinsic_curvature_gradient) =
          isotropic_prefactor * get<1>(x);
      get<2>(*trace_extrinsic_curvature_gradient) =
          isotropic_prefactor * get<2>(x);
      break;
    }
    case SchwarzschildCoordinates::MaximalIsotropic: {
      std::fill(trace_extrinsic_curvature_gradient->begin(),
                trace_extrinsic_curvature_gradient->end(), 0.);
      break;
    }
      // LCOV_EXCL_START
    default:
      ERROR("Missing case for SchwarzschildCoordinates");
      // LCOV_EXCL_END
  }
}

template <typename DataType>
void SchwarzschildVariables<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> conformal_factor_minus_one,
    const gsl::not_null<Cache*> cache,
    Xcts::Tags::ConformalFactorMinusOne<DataType> /*meta*/) const {
  switch (coordinate_system) {
    case SchwarzschildCoordinates::Isotropic: {
      const auto& r =
          get(cache->get_var(*this, detail::Tags::Radius<DataType>{}));
      get(*conformal_factor_minus_one) = 0.5 * mass / r;
      break;
    }
    case SchwarzschildCoordinates::PainleveGullstrand: {
      get(*conformal_factor_minus_one) = 0.;
      break;
    }
    case SchwarzschildCoordinates::KerrSchildIsotropic: {
      const auto& conformal_factor =
          get(cache->get_var(*this, Xcts::Tags::ConformalFactor<DataType>{}));
      get(*conformal_factor_minus_one) = conformal_factor - 1.;
      break;
    }
    case SchwarzschildCoordinates::MaximalIsotropic: {
      const auto& conformal_factor =
          get(cache->get_var(*this, Xcts::Tags::ConformalFactor<DataType>{}));
      get(*conformal_factor_minus_one) = conformal_factor - 1.;
      break;
    }
      // LCOV_EXCL_START
    default:
      ERROR("Missing case for SchwarzschildCoordinates");
      // LCOV_EXCL_END
  }
}

template <typename DataType>
void SchwarzschildVariables<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> conformal_factor,
    const gsl::not_null<Cache*> cache,
    Xcts::Tags::ConformalFactor<DataType> /*meta*/) const {
  switch (coordinate_system) {
    case SchwarzschildCoordinates::Isotropic: {
      const auto& conformal_factor_minus_one = get(cache->get_var(
          *this, Xcts::Tags::ConformalFactorMinusOne<DataType>{}));
      get(*conformal_factor) = conformal_factor_minus_one + 1.;
      break;
    }
    case SchwarzschildCoordinates::PainleveGullstrand: {
      get(*conformal_factor) = 1.;
      break;
    }
    case SchwarzschildCoordinates::KerrSchildIsotropic: {
      // Eq. (7.35) in https://arxiv.org/abs/gr-qc/0510016
      const auto& lapse =
          get(cache->get_var(*this, gr::Tags::Lapse<DataType>{}));
      get(*conformal_factor) = 2. * exp(1. / lapse - 1.) / (1. + 1. / lapse);
      break;
    }
    case SchwarzschildCoordinates::MaximalIsotropic: {
      const auto r_areal =
          get(cache->get_var(*this, detail::Tags::ArealRadius<DataType>{}));
      const auto& r_iso =
          get(cache->get_var(*this, detail::Tags::Radius<DataType>{}));
      get(*conformal_factor) = sqrt(r_areal / r_iso);
      break;
    }
      // LCOV_EXCL_START
    default:
      ERROR("Missing case for SchwarzschildCoordinates");
      // LCOV_EXCL_END
  }
}

template <typename DataType>
void SchwarzschildVariables<DataType>::operator()(
    const gsl::not_null<tnsr::i<DataType, 3>*> conformal_factor_gradient,
    const gsl::not_null<Cache*> cache,
    ::Tags::deriv<Xcts::Tags::ConformalFactorMinusOne<DataType>,
                  tmpl::size_t<3>, Frame::Inertial> /*meta*/) const {
  switch (coordinate_system) {
    case SchwarzschildCoordinates::Isotropic: {
      const auto& r =
          get(cache->get_var(*this, detail::Tags::Radius<DataType>{}));
      const DataType isotropic_prefactor = -0.5 * mass / cube(r);
      get<0>(*conformal_factor_gradient) = isotropic_prefactor * get<0>(x);
      get<1>(*conformal_factor_gradient) = isotropic_prefactor * get<1>(x);
      get<2>(*conformal_factor_gradient) = isotropic_prefactor * get<2>(x);
      break;
    }
    case SchwarzschildCoordinates::PainleveGullstrand: {
      std::fill(conformal_factor_gradient->begin(),
                conformal_factor_gradient->end(), 0.);
      break;
    }
    case SchwarzschildCoordinates::KerrSchildIsotropic: {
      const auto& rbar =
          get(cache->get_var(*this, detail::Tags::Radius<DataType>{}));
      const auto& r =
          get(cache->get_var(*this, detail::Tags::ArealRadius<DataType>{}));
      const auto& conformal_factor =
          get(cache->get_var(*this, Xcts::Tags::ConformalFactor<DataType>{}));
      const auto one_over_lapse = sqrt(1. + 2. * mass / r);
      const DataType isotropic_prefactor =
          -conformal_factor * mass /
          (r * (1. + one_over_lapse) * one_over_lapse * square(rbar));
      get<0>(*conformal_factor_gradient) = isotropic_prefactor * get<0>(x);
      get<1>(*conformal_factor_gradient) = isotropic_prefactor * get<1>(x);
      get<2>(*conformal_factor_gradient) = isotropic_prefactor * get<2>(x);
      break;
    }
    case SchwarzschildCoordinates::MaximalIsotropic: {
      const auto& r_iso =
          get(cache->get_var(*this, detail::Tags::Radius<DataType>{}));
      const auto& r_areal =
          get(cache->get_var(*this, detail::Tags::ArealRadius<DataType>{}));

      const DataType drdR =
          maximal_isotropic_radius_from_areal_deriv(r_areal, r_iso, mass);
      const DataType dPsidR = .5 * sqrt(r_iso / r_areal) *
                              (1. / r_iso - (r_areal / square(r_iso)) * drdR);
      const DataType prefactor = dPsidR / (drdR * r_iso);
      get<0>(*conformal_factor_gradient) = prefactor * get<0>(x);
      get<1>(*conformal_factor_gradient) = prefactor * get<1>(x);
      get<2>(*conformal_factor_gradient) = prefactor * get<2>(x);
      break;
    }
      // LCOV_EXCL_START
    default:
      ERROR("Missing case for SchwarzschildCoordinates");
      // LCOV_EXCL_END
  }
}

template <typename DataType>
void SchwarzschildVariables<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*>
        lapse_times_conformal_factor_minus_one,
    const gsl::not_null<Cache*> cache,
    Xcts::Tags::LapseTimesConformalFactorMinusOne<DataType> /*meta*/) const {
  switch (coordinate_system) {
    case SchwarzschildCoordinates::Isotropic: {
      const auto& r =
          get(cache->get_var(*this, detail::Tags::Radius<DataType>{}));
      get(*lapse_times_conformal_factor_minus_one) = -0.5 * mass / r;
      break;
    }
    case SchwarzschildCoordinates::PainleveGullstrand: {
      get(*lapse_times_conformal_factor_minus_one) = 0.;
      break;
    }
    case SchwarzschildCoordinates::KerrSchildIsotropic: {
      const auto& lapse_times_conformal_factor = get(cache->get_var(
          *this, Xcts::Tags::LapseTimesConformalFactor<DataType>{}));
      get(*lapse_times_conformal_factor_minus_one) =
          lapse_times_conformal_factor - 1.;
      break;
    }
    case SchwarzschildCoordinates::MaximalIsotropic: {
      const auto& lapse_times_conformal_factor = get(cache->get_var(
          *this, Xcts::Tags::LapseTimesConformalFactor<DataType>{}));
      get(*lapse_times_conformal_factor_minus_one) =
          lapse_times_conformal_factor - 1.;
      break;
    }
      // LCOV_EXCL_START
    default:
      ERROR("Missing case for SchwarzschildCoordinates");
      // LCOV_EXCL_END
  }
}

template <typename DataType>
void SchwarzschildVariables<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> lapse_times_conformal_factor,
    const gsl::not_null<Cache*> cache,
    Xcts::Tags::LapseTimesConformalFactor<DataType> /*meta*/) const {
  switch (coordinate_system) {
    case SchwarzschildCoordinates::Isotropic: {
      const auto& lapse_times_conformal_factor_minus_one = get(cache->get_var(
          *this, Xcts::Tags::LapseTimesConformalFactorMinusOne<DataType>{}));
      get(*lapse_times_conformal_factor) =
          lapse_times_conformal_factor_minus_one + 1.;
      break;
    }
    case SchwarzschildCoordinates::PainleveGullstrand: {
      get(*lapse_times_conformal_factor) = 1.;
      break;
    }
    case SchwarzschildCoordinates::KerrSchildIsotropic: {
      const auto& lapse =
          get(cache->get_var(*this, gr::Tags::Lapse<DataType>{}));
      const auto& conformal_factor =
          get(cache->get_var(*this, Xcts::Tags::ConformalFactor<DataType>{}));
      get(*lapse_times_conformal_factor) = lapse * conformal_factor;
      break;
    }
    case SchwarzschildCoordinates::MaximalIsotropic: {
      const auto& lapse =
          get(cache->get_var(*this, gr::Tags::Lapse<DataType>{}));
      const auto& conformal_factor =
          get(cache->get_var(*this, Xcts::Tags::ConformalFactor<DataType>{}));
      get(*lapse_times_conformal_factor) = lapse * conformal_factor;
      break;
    }
      // LCOV_EXCL_START
    default:
      ERROR("Missing case for SchwarzschildCoordinates");
      // LCOV_EXCL_END
  }
}

template <typename DataType>
void SchwarzschildVariables<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> lapse,
    const gsl::not_null<Cache*> cache,
    gr::Tags::Lapse<DataType> /*meta*/) const {
  switch (coordinate_system) {
    case SchwarzschildCoordinates::Isotropic: {
      const auto& conformal_factor =
          get(cache->get_var(*this, Xcts::Tags::ConformalFactor<DataType>{}));
      const auto& lapse_times_conformal_factor = get(cache->get_var(
          *this, Xcts::Tags::LapseTimesConformalFactor<DataType>{}));
      get(*lapse) = lapse_times_conformal_factor / conformal_factor;
      break;
    }
    case SchwarzschildCoordinates::PainleveGullstrand: {
      get(*lapse) = 1.;
      break;
    }
    case SchwarzschildCoordinates::KerrSchildIsotropic: {
      const auto& r =
          get(cache->get_var(*this, detail::Tags::ArealRadius<DataType>{}));
      get(*lapse) = 1. / sqrt(1. + 2. * mass / r);
      break;
    }
    case SchwarzschildCoordinates::MaximalIsotropic: {
      const auto& r_areal =
          get(cache->get_var(*this, detail::Tags::ArealRadius<DataType>{}));
      get(*lapse) =
          sqrt(1. - 2. * mass / r_areal +
               27. * square(square(mass)) / (16. * square(square(r_areal))));
      break;
    }
      // LCOV_EXCL_START
    default:
      ERROR("Missing case for SchwarzschildCoordinates");
      // LCOV_EXCL_END
  }
}

template <typename DataType>
void SchwarzschildVariables<DataType>::operator()(
    const gsl::not_null<tnsr::i<DataType, 3>*> deriv_lapse,
    const gsl::not_null<Cache*> cache,
    ::Tags::deriv<gr::Tags::Lapse<DataType>, tmpl::size_t<3>,
                  Frame::Inertial> /*meta*/) const {
  switch (coordinate_system) {
    case SchwarzschildCoordinates::Isotropic: {
      const auto& r =
          get(cache->get_var(*this, detail::Tags::Radius<DataType>{}));
      const auto& conformal_factor =
          cache->get_var(*this, Xcts::Tags::ConformalFactor<DataType>{});
      const DataType isotropic_prefactor =
          mass / square(get(conformal_factor)) / cube(r);
      get<0>(*deriv_lapse) = isotropic_prefactor * get<0>(x);
      get<1>(*deriv_lapse) = isotropic_prefactor * get<1>(x);
      get<2>(*deriv_lapse) = isotropic_prefactor * get<2>(x);
      break;
    }
    case SchwarzschildCoordinates::PainleveGullstrand: {
      std::fill(deriv_lapse->begin(), deriv_lapse->end(), 0.);
      break;
    }
    case SchwarzschildCoordinates::KerrSchildIsotropic: {
      const auto& r =
          get(cache->get_var(*this, detail::Tags::ArealRadius<DataType>{}));
      const auto& rbar =
          get(cache->get_var(*this, detail::Tags::Radius<DataType>{}));
      const auto& lapse =
          get(cache->get_var(*this, gr::Tags::Lapse<DataType>{}));
      const DataType isotropic_prefactor =
          mass * pow<4>(lapse) / (r * square(rbar));
      get<0>(*deriv_lapse) = isotropic_prefactor * get<0>(x);
      get<1>(*deriv_lapse) = isotropic_prefactor * get<1>(x);
      get<2>(*deriv_lapse) = isotropic_prefactor * get<2>(x);
      break;
    }
    case SchwarzschildCoordinates::MaximalIsotropic: {
      const auto& r_areal =
          get(cache->get_var(*this, detail::Tags::ArealRadius<DataType>{}));
      const auto& r_iso =
          get(cache->get_var(*this, detail::Tags::Radius<DataType>{}));
      const auto& lapse =
          get(cache->get_var(*this, gr::Tags::Lapse<DataType>{}));

      const DataType drdR =
          maximal_isotropic_radius_from_areal_deriv(r_areal, r_iso, mass);
      const DataType dLapsedR =
          (1. / (2. * lapse)) *
          (2. * mass / square(r_areal) -
           27. * square(square(mass)) / (4. * pow(r_areal, 5)));
      const DataType prefactor = dLapsedR / (drdR * r_iso);
      get<0>(*deriv_lapse) = prefactor * get<0>(x);
      get<1>(*deriv_lapse) = prefactor * get<1>(x);
      get<2>(*deriv_lapse) = prefactor * get<2>(x);
      break;
    }
      // LCOV_EXCL_START
    default:
      ERROR("Missing case for SchwarzschildCoordinates");
      // LCOV_EXCL_END
  }
}

template <typename DataType>
void SchwarzschildVariables<DataType>::operator()(
    const gsl::not_null<tnsr::i<DataType, 3>*>
        lapse_times_conformal_factor_gradient,
    const gsl::not_null<Cache*> cache,
    ::Tags::deriv<Xcts::Tags::LapseTimesConformalFactorMinusOne<DataType>,
                  tmpl::size_t<3>, Frame::Inertial> /*meta*/) const {
  switch (coordinate_system) {
    case SchwarzschildCoordinates::Isotropic: {
      *lapse_times_conformal_factor_gradient = cache->get_var(
          *this, ::Tags::deriv<Xcts::Tags::ConformalFactorMinusOne<DataType>,
                               tmpl::size_t<3>, Frame::Inertial>{});
      get<0>(*lapse_times_conformal_factor_gradient) *= -1.;
      get<1>(*lapse_times_conformal_factor_gradient) *= -1.;
      get<2>(*lapse_times_conformal_factor_gradient) *= -1.;
      break;
    }
    case SchwarzschildCoordinates::PainleveGullstrand: {
      std::fill(lapse_times_conformal_factor_gradient->begin(),
                lapse_times_conformal_factor_gradient->end(), 0.);
      break;
    }
    case SchwarzschildCoordinates::KerrSchildIsotropic: {
      const auto& rbar =
          get(cache->get_var(*this, detail::Tags::Radius<DataType>{}));
      const auto& r =
          get(cache->get_var(*this, detail::Tags::ArealRadius<DataType>{}));
      const auto& conformal_factor =
          get(cache->get_var(*this, Xcts::Tags::ConformalFactor<DataType>{}));
      const auto& conformal_factor_gradient = cache->get_var(
          *this, ::Tags::deriv<Xcts::Tags::ConformalFactorMinusOne<DataType>,
                               tmpl::size_t<3>, Frame::Inertial>{});
      const auto& lapse =
          get(cache->get_var(*this, gr::Tags::Lapse<DataType>{}));
      const DataType isotropic_prefactor =
          conformal_factor * pow<4>(lapse) * mass / r / square(rbar);
      *lapse_times_conformal_factor_gradient = conformal_factor_gradient;
      get<0>(*lapse_times_conformal_factor_gradient) *= lapse;
      get<1>(*lapse_times_conformal_factor_gradient) *= lapse;
      get<2>(*lapse_times_conformal_factor_gradient) *= lapse;
      get<0>(*lapse_times_conformal_factor_gradient) +=
          isotropic_prefactor * get<0>(x);
      get<1>(*lapse_times_conformal_factor_gradient) +=
          isotropic_prefactor * get<1>(x);
      get<2>(*lapse_times_conformal_factor_gradient) +=
          isotropic_prefactor * get<2>(x);
      break;
    }
    case SchwarzschildCoordinates::MaximalIsotropic: {
      const auto& conformal_factor =
          get(cache->get_var(*this, Xcts::Tags::ConformalFactor<DataType>{}));
      const auto& conformal_factor_gradient = cache->get_var(
          *this, ::Tags::deriv<Xcts::Tags::ConformalFactorMinusOne<DataType>,
                               tmpl::size_t<3>, Frame::Inertial>{});
      const auto& lapse =
          get(cache->get_var(*this, gr::Tags::Lapse<DataType>{}));
      const auto& deriv_lapse = cache->get_var(
          *this, ::Tags::deriv<gr::Tags::Lapse<DataType>, tmpl::size_t<3>,
                               Frame::Inertial>{});
      *lapse_times_conformal_factor_gradient = conformal_factor_gradient;
      get<0>(*lapse_times_conformal_factor_gradient) *= lapse;
      get<1>(*lapse_times_conformal_factor_gradient) *= lapse;
      get<2>(*lapse_times_conformal_factor_gradient) *= lapse;
      get<0>(*lapse_times_conformal_factor_gradient) +=
          conformal_factor * get<0>(deriv_lapse);
      get<1>(*lapse_times_conformal_factor_gradient) +=
          conformal_factor * get<1>(deriv_lapse);
      get<2>(*lapse_times_conformal_factor_gradient) +=
          conformal_factor * get<2>(deriv_lapse);
      break;
    }
      // LCOV_EXCL_START
    default:
      ERROR("Missing case for SchwarzschildCoordinates");
      // LCOV_EXCL_END
  }
}

// Set the background shift to zero in the decomposition:
//   shift = shift_background + shift_excess
// See docs of Xcts::Tags::ShiftExcess.
template <typename DataType>
void SchwarzschildVariables<DataType>::operator()(
    const gsl::not_null<tnsr::I<DataType, 3>*> shift_background,
    const gsl::not_null<Cache*> /*cache*/,
    Xcts::Tags::ShiftBackground<DataType, 3, Frame::Inertial> /*meta*/) const {
  std::fill(shift_background->begin(), shift_background->end(), 0.);
}

template <typename DataType>
void SchwarzschildVariables<DataType>::operator()(
    const gsl::not_null<tnsr::II<DataType, 3, Frame::Inertial>*>
        longitudinal_shift_background_minus_dt_conformal_metric,
    const gsl::not_null<Cache*> /*cache*/,
    Xcts::Tags::LongitudinalShiftBackgroundMinusDtConformalMetric<
        DataType, 3, Frame::Inertial> /*meta*/) const {
  std::fill(longitudinal_shift_background_minus_dt_conformal_metric->begin(),
            longitudinal_shift_background_minus_dt_conformal_metric->end(), 0.);
}

template <typename DataType>
void SchwarzschildVariables<DataType>::operator()(
    const gsl::not_null<tnsr::I<DataType, 3>*> shift_excess,
    const gsl::not_null<Cache*> cache,
    Xcts::Tags::ShiftExcess<DataType, 3, Frame::Inertial> /*meta*/) const {
  switch (coordinate_system) {
    case SchwarzschildCoordinates::Isotropic: {
      std::fill(shift_excess->begin(), shift_excess->end(), 0.);
      break;
    }
    case SchwarzschildCoordinates::PainleveGullstrand: {
      const auto& r =
          get(cache->get_var(*this, detail::Tags::Radius<DataType>{}));
      const DataType isotropic_prefactor = sqrt(2. * mass) / pow(r, 1.5);
      *shift_excess = x;
      get<0>(*shift_excess) *= isotropic_prefactor;
      get<1>(*shift_excess) *= isotropic_prefactor;
      get<2>(*shift_excess) *= isotropic_prefactor;
      break;
    }
    case SchwarzschildCoordinates::KerrSchildIsotropic: {
      const auto& r =
          get(cache->get_var(*this, detail::Tags::ArealRadius<DataType>{}));
      const auto& lapse =
          get(cache->get_var(*this, gr::Tags::Lapse<DataType>{}));
      const DataType isotropic_prefactor = 2. * mass * lapse / square(r);
      *shift_excess = x;
      get<0>(*shift_excess) *= isotropic_prefactor;
      get<1>(*shift_excess) *= isotropic_prefactor;
      get<2>(*shift_excess) *= isotropic_prefactor;
      break;
    }
    case SchwarzschildCoordinates::MaximalIsotropic: {
      const auto& r_areal =
          get(cache->get_var(*this, detail::Tags::ArealRadius<DataType>{}));
      get<0>(*shift_excess) =
          .75 * sqrt(3) * square(mass) * get<0>(x) / pow(r_areal, 3);
      get<1>(*shift_excess) =
          .75 * sqrt(3) * square(mass) * get<1>(x) / pow(r_areal, 3);
      get<2>(*shift_excess) =
          .75 * sqrt(3) * square(mass) * get<2>(x) / pow(r_areal, 3);
      break;
    }
      // LCOV_EXCL_START
    default:
      ERROR("Missing case for SchwarzschildCoordinates");
      // LCOV_EXCL_END
  }
}

template <typename DataType>
void SchwarzschildVariables<DataType>::operator()(
    const gsl::not_null<tnsr::iJ<DataType, 3>*> deriv_shift_excess,
    const gsl::not_null<Cache*> cache,
    ::Tags::deriv<Xcts::Tags::ShiftExcess<DataType, 3, Frame::Inertial>,
                  tmpl::size_t<3>, Frame::Inertial> /*meta*/) const {
  switch (coordinate_system) {
    case SchwarzschildCoordinates::Isotropic: {
      std::fill(deriv_shift_excess->begin(), deriv_shift_excess->end(), 0.);
      break;
    }
    case SchwarzschildCoordinates::PainleveGullstrand: {
      const auto& r =
          get(cache->get_var(*this, detail::Tags::Radius<DataType>{}));
      const DataType diagonal_prefactor = sqrt(2. * mass) / pow(r, 1.5);
      const DataType isotropic_prefactor =
          -1.5 * diagonal_prefactor / square(r);
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          deriv_shift_excess->get(i, j) =
              isotropic_prefactor * x.get(i) * x.get(j);
        }
        deriv_shift_excess->get(i, i) += diagonal_prefactor;
      }
      break;
    }
    case SchwarzschildCoordinates::KerrSchildIsotropic: {
      const auto& rbar =
          get(cache->get_var(*this, detail::Tags::Radius<DataType>{}));
      const auto& r =
          get(cache->get_var(*this, detail::Tags::ArealRadius<DataType>{}));
      const auto& lapse =
          get(cache->get_var(*this, gr::Tags::Lapse<DataType>{}));
      const DataType diagonal_prefactor = 2. * mass * lapse / square(r);
      const DataType isotropic_prefactor = diagonal_prefactor *
                                           (square(lapse) * mass / r - 2.) *
                                           lapse / square(rbar);
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          deriv_shift_excess->get(i, j) =
              isotropic_prefactor * x.get(i) * x.get(j);
        }
        deriv_shift_excess->get(i, i) += diagonal_prefactor;
      }
      break;
    }
    case SchwarzschildCoordinates::MaximalIsotropic: {
      const auto& r_areal =
          get(cache->get_var(*this, detail::Tags::ArealRadius<DataType>{}));
      const auto& r_iso =
          get(cache->get_var(*this, detail::Tags::Radius<DataType>{}));
      const DataType BetaRadial =
          .75 * sqrt(3) * square(mass) * r_iso / pow(r_areal, 3);

      const DataType drdR =
          maximal_isotropic_radius_from_areal_deriv(r_areal, r_iso, mass);
      const DataType dBetaRadialdR =
          .75 * sqrt(3) * square(mass) *
          (drdR / pow(r_areal, 3) - 3. * r_iso / pow(r_areal, 4));
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j < 3; ++j) {
          deriv_shift_excess->get(i, j) =
              (dBetaRadialdR / (drdR * square(r_iso)) -
               BetaRadial / pow(r_iso, 3)) *
              x.get(i) * x.get(j);
        }
        deriv_shift_excess->get(i, i) += BetaRadial / r_iso;
      }
      break;
    }
      // LCOV_EXCL_START
    default:
      ERROR("Missing case for SchwarzschildCoordinates");
      // LCOV_EXCL_END
  }
}

template <typename DataType>
void SchwarzschildVariables<DataType>::operator()(
    const gsl::not_null<tnsr::ii<DataType, 3>*> extrinsic_curvature,
    const gsl::not_null<Cache*> cache,
    gr::Tags::ExtrinsicCurvature<DataType, 3> /*meta*/) const {
  switch (coordinate_system) {
    case SchwarzschildCoordinates::Isotropic: {
      std::fill(extrinsic_curvature->begin(), extrinsic_curvature->end(), 0.);
      break;
    }
    case SchwarzschildCoordinates::PainleveGullstrand: {
      const auto& r =
          get(cache->get_var(*this, detail::Tags::Radius<DataType>{}));
      const DataType diagonal_prefactor = sqrt(2. * mass / cube(r));
      const DataType isotropic_prefactor =
          -3. / 2. * diagonal_prefactor / square(r);
      for (size_t i = 0; i < 3; ++i) {
        for (size_t j = 0; j <= i; ++j) {
          extrinsic_curvature->get(i, j) =
              isotropic_prefactor * x.get(i) * x.get(j);
        }
        extrinsic_curvature->get(i, i) += diagonal_prefactor;
      }
      break;
    }
    case SchwarzschildCoordinates::KerrSchildIsotropic: {
      // Background shift and \bar{u} are both zero
      const auto& longitudinal_shift_minus_dt_conformal_metric = cache->get_var(
          *this,
          Xcts::Tags::LongitudinalShiftExcess<DataType, 3, Frame::Inertial>{});
      Xcts::extrinsic_curvature(
          extrinsic_curvature,
          cache->get_var(*this, Xcts::Tags::ConformalFactor<DataType>{}),
          cache->get_var(*this, gr::Tags::Lapse<DataType>{}),
          cache->get_var(
              *this,
              Xcts::Tags::ConformalMetric<DataType, 3, Frame::Inertial>{}),
          longitudinal_shift_minus_dt_conformal_metric,
          cache->get_var(*this, gr::Tags::TraceExtrinsicCurvature<DataType>{}));
      break;
    }
    case SchwarzschildCoordinates::MaximalIsotropic: {
      // Background shift and \bar{u} are both zero
      const auto& longitudinal_shift_minus_dt_conformal_metric = cache->get_var(
          *this,
          Xcts::Tags::LongitudinalShiftExcess<DataType, 3, Frame::Inertial>{});
      Xcts::extrinsic_curvature(
          extrinsic_curvature,
          cache->get_var(*this, Xcts::Tags::ConformalFactor<DataType>{}),
          cache->get_var(*this, gr::Tags::Lapse<DataType>{}),
          cache->get_var(
              *this,
              Xcts::Tags::ConformalMetric<DataType, 3, Frame::Inertial>{}),
          longitudinal_shift_minus_dt_conformal_metric,
          cache->get_var(*this, gr::Tags::TraceExtrinsicCurvature<DataType>{}));
      break;
    }
      // LCOV_EXCL_START
    default:
      ERROR("Missing case for SchwarzschildCoordinates");
      // LCOV_EXCL_END
  }
}

template <typename DataType>
void SchwarzschildVariables<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> energy_density,
    const gsl::not_null<Cache*> /*cache*/,
    gr::Tags::Conformal<gr::Tags::EnergyDensity<DataType>,
                        ConformalMatterScale> /*meta*/) const {
  std::fill(energy_density->begin(), energy_density->end(), 0.);
}

template <typename DataType>
void SchwarzschildVariables<DataType>::operator()(
    const gsl::not_null<Scalar<DataType>*> stress_trace,
    const gsl::not_null<Cache*> /*cache*/,
    gr::Tags::Conformal<gr::Tags::StressTrace<DataType>,
                        ConformalMatterScale> /*meta*/) const {
  std::fill(stress_trace->begin(), stress_trace->end(), 0.);
}

template <typename DataType>
void SchwarzschildVariables<DataType>::operator()(
    const gsl::not_null<tnsr::I<DataType, 3>*> momentum_density,
    const gsl::not_null<Cache*> /*cache*/,
    gr::Tags::Conformal<gr::Tags::MomentumDensity<DataType, 3>,
                        ConformalMatterScale> /*meta*/) const {
  std::fill(momentum_density->begin(), momentum_density->end(), 0.);
}

template class SchwarzschildVariables<double>;
template class SchwarzschildVariables<DataVector>;

}  // namespace detail

PUP::able::PUP_ID Schwarzschild::my_PUP_ID = 0;  // NOLINT

}  // namespace Xcts::Solutions

// Instantiate implementations for common variables
template class Xcts::Solutions::CommonVariables<
    double,
    typename Xcts::Solutions::detail::SchwarzschildVariables<double>::Cache>;
template class Xcts::Solutions::CommonVariables<
    DataVector, typename Xcts::Solutions::detail::SchwarzschildVariables<
                    DataVector>::Cache>;
template class Xcts::AnalyticData::CommonVariables<
    double,
    typename Xcts::Solutions::detail::SchwarzschildVariables<double>::Cache>;
template class Xcts::AnalyticData::CommonVariables<
    DataVector, typename Xcts::Solutions::detail::SchwarzschildVariables<
                    DataVector>::Cache>;
