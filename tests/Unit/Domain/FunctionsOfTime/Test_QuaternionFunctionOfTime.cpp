// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <catch2/matchers/catch_matchers.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/FunctionsOfTime/QuaternionFunctionOfTime.hpp"
#include "Framework/TestHelpers.hpp"
#include "Helpers/Utilities/Serialization/Versioning.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/Serialization/Serialize.hpp"

namespace {
void test_serialization_versioning() {
  using QuatFoT = domain::FunctionsOfTime::QuaternionFunctionOfTime<2>;
  register_classes_with_charm<QuatFoT>();
  const std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime> func(
      std::make_unique<QuatFoT>(
          0.0, std::array<DataVector, 1>{DataVector{{1.0, 0.0, 0.0, 0.0}}},
          std::array<DataVector, 3>{DataVector{3, 0.0},
                                    DataVector{0.0, 0.0, 3.78},
                                    DataVector{3, 0.0}},
          0.6));
  func->update(0.6, DataVector{3, 0.0}, 1.0);

  TestHelpers::serialization::test_versioning<QuatFoT>(
      "Domain/FunctionsOfTime/QuaternionFunctionOfTime.serializations",
      "version 5", func);
}

void test_out_of_order_update() {
  using QuatFoT = domain::FunctionsOfTime::QuaternionFunctionOfTime<2>;
  register_classes_with_charm<QuatFoT>();
  const DataVector quat_dv{1.0, 0.0, 0.0, 0.0};
  const DataVector zero_dv{3, 0.0};
  const std::unique_ptr<domain::FunctionsOfTime::FunctionOfTime> func(
      std::make_unique<QuatFoT>(
          0.0, std::array<DataVector, 1>{quat_dv},
          std::array<DataVector, 3>{zero_dv, zero_dv, zero_dv}, 0.5));

  func->update(2.5, zero_dv, 3.0);
  func->update(1.5, zero_dv, 2.0);

  CHECK_THROWS_WITH(
      func->func(0.7),
      Catch::Matchers::ContainsSubstring("Attempt to evaluate at time") and
          Catch::Matchers::ContainsSubstring(
              ", which is after the expiration time 5"));
  CHECK_THROWS_WITH(
      func->func(1.3),
      Catch::Matchers::ContainsSubstring("Attempt to evaluate at time 1.30") and
          Catch::Matchers::ContainsSubstring(
              ", which is after the expiration time 5"));
  CHECK_THROWS_WITH(
      func->func(1.9),
      Catch::Matchers::ContainsSubstring("Attempt to evaluate at time") and
          Catch::Matchers::ContainsSubstring(
              ", which is after the expiration time 5"));

  func->update(0.5, zero_dv, 1.0);

  // This one should pass now
  CHECK(func->func(0.7) == std::array{quat_dv});
  CHECK_THROWS_WITH(
      func->func(1.3),
      Catch::Matchers::ContainsSubstring("Attempt to evaluate at time 1.30") and
          Catch::Matchers::ContainsSubstring(
              ", which is after the expiration time 1.0"));
  CHECK_THROWS_WITH(
      func->func(1.9),
      Catch::Matchers::ContainsSubstring("Attempt to evaluate at time") and
          Catch::Matchers::ContainsSubstring(
              ", which is after the expiration time 1.0"));

  func->update(2.0, zero_dv, 2.5);

  CHECK_THROWS_WITH(
      func->func(1.3),
      Catch::Matchers::ContainsSubstring("Attempt to evaluate at time 1.30") and
          Catch::Matchers::ContainsSubstring(
              ", which is after the expiration time 1.0"));
  CHECK_THROWS_WITH(
      func->func(1.9),
      Catch::Matchers::ContainsSubstring("Attempt to evaluate at time") and
          Catch::Matchers::ContainsSubstring(
              ", which is after the expiration time 1.0"));

  func->update(1.0, zero_dv, 1.5);

  CHECK(func->func(1.3) == std::array{quat_dv});
  CHECK(func->func(1.9) == std::array{quat_dv});
  CHECK(func->func(2.2) == std::array{quat_dv});
  CHECK(func->func(2.8) == std::array{quat_dv});
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Domain.FunctionsOfTime.QuaternionFunctionOfTime",
                  "[Unit][Domain]") {
  {
    INFO("QuaternionFunctionOfTime: Time bounds");
    domain::FunctionsOfTime::QuaternionFunctionOfTime<2> qfot{
        0.0, std::array<DataVector, 1>{DataVector{4, 0.0}},
        std::array<DataVector, 3>{DataVector{3, 0.0}, DataVector{3, 0.0},
                                  DataVector{3, 0.0}},
        0.5};
    CHECK(qfot.time_bounds() == std::array<double, 2>({0.0, 0.5}));
  }

  {
    INFO("QuaternionFunctionOfTime: Check output");
    domain::FunctionsOfTime::QuaternionFunctionOfTime<2> qfot{
        -0.1, std::array<DataVector, 1>{DataVector{{1.0, 0.0, 0.0, 0.0}}},
        std::array<DataVector, 3>{DataVector{{0.0, 0.0, -0.3}},
                                  DataVector{0.0, 0.0, 0.145},
                                  DataVector{3, 0.0}},
        0.5};

    CHECK(get_output(domain::FunctionsOfTime::QuaternionFunctionOfTime<2>{}) ==
          "Quaternion:\n"
          "backlog=()\n"
          "Angle:\n"
          "backlog=()");
    const std::string expected_output =
        "Quaternion:\n"
        "t=-0.1: (1,0,0,0)\n"
        "backlog=()\n"
        "Angle:\n"
        "t=-0.1: (0,0,-0.3) (0,0,0.145) (0,0,0)\n"
        "backlog=()";

    const std::string output = get_output(qfot);
    CHECK(output == expected_output);
  }

  {
    INFO("QuaternionFunctionOfTime: Internal PiecewisePolynomial");
    DataVector init_omega{0.0, 0.0, 3.78};
    domain::FunctionsOfTime::QuaternionFunctionOfTime<2> qfot{
        0.0, std::array<DataVector, 1>{DataVector{{1.0, 0.0, 0.0, 0.0}}},
        std::array<DataVector, 3>{DataVector{3, 0.0}, init_omega,
                                  DataVector{3, 0.0}},
        0.6};
    domain::FunctionsOfTime::PiecewisePolynomial<2> pp{
        0.0,
        std::array<DataVector, 3>{DataVector{3, 0.0}, init_omega,
                                  DataVector{3, 0.0}},
        0.6};
    qfot.update(0.6, DataVector{3, 0.0}, 1.0);
    pp.update(0.6, DataVector{3, 0.0}, 1.0);

    CHECK(qfot.angle_func(0.4) == pp.func(0.4));
    CHECK(qfot.angle_func_and_deriv(0.4) == pp.func_and_deriv(0.4));
    CHECK(qfot.angle_func_and_2_derivs(0.4) == pp.func_and_2_derivs(0.4));
  }

  {
    INFO("QuaternionFunctionOfTime: pup, cloning, extra functions, copy/move");
    DataVector init_omega{0.0, 0.0, 1.0};
    domain::FunctionsOfTime::QuaternionFunctionOfTime<2> qfot{
        0.0, std::array<DataVector, 1>{DataVector{{1.0, 0.0, 0.0, 0.0}}},
        std::array<DataVector, 3>{DataVector{3, 0.0}, init_omega,
                                  DataVector{3, 0.0}},
        2.5};
    // Different expiration time to check comparison operators
    domain::FunctionsOfTime::QuaternionFunctionOfTime<2> qfot2{
        0.0, std::array<DataVector, 1>{DataVector{{1.0, 0.0, 0.0, 0.0}}},
        std::array<DataVector, 3>{DataVector{3, 0.0}, init_omega,
                                  DataVector{3, 0.0}},
        3.0};

    auto qfot_ptr = qfot.get_clone();
    domain::FunctionsOfTime::QuaternionFunctionOfTime<2>
        qfot_serialized_deserialized = serialize_and_deserialize(qfot);

    CHECK(qfot == qfot_serialized_deserialized);
    CHECK(qfot != qfot2);

    std::array<DataVector, 1> expected_func{
        DataVector{{cos(1.0), 0.0, 0.0, sin(1.0)}}};
    std::array<DataVector, 2> expected_func_and_deriv{
        DataVector{{cos(1.0), 0.0, 0.0, sin(1.0)}},
        DataVector{{-0.5 * sin(1.0), 0.0, 0.0, 0.5 * cos(1.0)}}};

    const auto do_check = [&](const auto& qfot_to_check,
                              const auto& ptr_to_check,
                              const auto& serialized_to_check) {
      CHECK_ITERABLE_APPROX(qfot_to_check.quat_func(2.0), expected_func);
      CHECK_ITERABLE_APPROX(qfot_to_check.quat_func_and_deriv(2.0),
                            expected_func_and_deriv);
      CHECK_ITERABLE_APPROX(qfot_to_check.func(2.0), expected_func);
      CHECK_ITERABLE_APPROX(qfot_to_check.func_and_deriv(2.0),
                            expected_func_and_deriv);
      CHECK_ITERABLE_APPROX(ptr_to_check->func(2.0), qfot_to_check.func(2.0));
      CHECK_ITERABLE_APPROX(ptr_to_check->func_and_deriv(2.0),
                            qfot_to_check.func_and_deriv(2.0));
      CHECK_ITERABLE_APPROX(serialized_to_check.func(2.0),
                            qfot_to_check.func(2.0));
      CHECK_ITERABLE_APPROX(serialized_to_check.func_and_deriv(2.0),
                            qfot_to_check.func_and_deriv(2.0));
    };

    do_check(qfot, qfot_ptr, qfot_serialized_deserialized);

    CHECK_THROWS_WITH(
        qfot.quat_func_and_3_derivs(2.0),
        Catch::Matchers::ContainsSubstring(
            "Need more angle derivs to compute the third derivative of the "
            "quaternion. Currently only have 2"));

    test_copy_semantics(qfot);
    test_move_semantics(std::move(qfot_serialized_deserialized), qfot);

    const auto copy_at_time = qfot.create_at_time(0.7, 2.8);
    const auto clone_copy_at_time = copy_at_time->get_clone();
    const auto copy_serialized = serialize_and_deserialize(
        dynamic_cast<
            const domain::FunctionsOfTime::QuaternionFunctionOfTime<2>&>(
            *copy_at_time));

    do_check(dynamic_cast<
                 const domain::FunctionsOfTime::QuaternionFunctionOfTime<2>&>(
                 *copy_at_time),
             clone_copy_at_time, copy_serialized);
  }

  {
    INFO("QuaternionFunctionOfTime: Constant omega");
    double t = 0.0;
    double expir_time = 0.5;
    const double omega_z = 1.3333;
    DataVector init_quat{1.0, 0.0, 0.0, 0.0};
    DataVector init_angle{0.0, 0.0, 0.0};
    DataVector init_omega{0.0, 0.0, omega_z};
    DataVector init_dtomega{0.0, 0.0, 0.0};

    // Construct QuaternionFunctionOfTime
    domain::FunctionsOfTime::QuaternionFunctionOfTime<2> qfot{
        t, std::array<DataVector, 1>{init_quat},
        std::array<DataVector, 3>{init_angle, init_omega, init_dtomega},
        expir_time};

    // Update stored PiecewisePolynomial with 0 2nd derivative so it's
    // constant. This will automatically update the stored quaternions as well
    const double time_step = 0.5;
    for (int i = 0; i < 15; i++) {
      t += time_step;
      expir_time += time_step;
      qfot.update(t, DataVector{3, 0.0}, expir_time);
      CHECK(qfot.expiration_after(t) == expir_time);
      CHECK(qfot.expiration_after(t + 0.5 * time_step) == expir_time);
      CHECK(qfot.expiration_after(t - 0.5 * time_step) == t);
    }
    CHECK(qfot.expiration_after(0.0) == time_step);

    // Get the quaternion and 2 derivatives at a certain time.
    double check_time = 5.398;
    const std::array<DataVector, 3> quat_func_and_2_derivs =
        qfot.quat_func_and_2_derivs(check_time);
    const std::array<DataVector, 3> quat_func_and_2_derivs2 =
        qfot.func_and_2_derivs(check_time);
    for (size_t i = 0; i < 3; i++) {
      CHECK_ITERABLE_APPROX(gsl::at(quat_func_and_2_derivs, i),
                            gsl::at(quat_func_and_2_derivs2, i));
    }

    // Analytic solution for constant omega
    // quat = ( cos(omega*t/2), 0, 0, sin(omega*t/2) )
    DataVector a_quat{{cos(0.5 * omega_z * check_time), 0.0, 0.0,
                       sin(0.5 * omega_z * check_time)}};
    DataVector a_dtquat{{-0.5 * omega_z * sin(0.5 * omega_z * check_time), 0.0,
                         0.0, 0.5 * omega_z * cos(0.5 * omega_z * check_time)}};
    DataVector a_dt2quat{
        {-0.25 * omega_z * omega_z * cos(0.5 * omega_z * check_time), 0.0, 0.0,
         -0.25 * omega_z * omega_z * sin(0.5 * omega_z * check_time)}};

    // Compare analytic solution to numerical
    Approx custom_approx = Approx::custom().epsilon(1.0e-12).scale(1.0);
    {
      INFO("  Compare quaternion");
      CHECK_ITERABLE_CUSTOM_APPROX(quat_func_and_2_derivs[0], a_quat,
                                   custom_approx);
    }
    {
      INFO("  Compare derivative of quaternion");
      CHECK_ITERABLE_CUSTOM_APPROX(quat_func_and_2_derivs[1], a_dtquat,
                                   custom_approx);
    }
    {
      INFO("  Compare second derivative of quaternion");
      CHECK_ITERABLE_CUSTOM_APPROX(quat_func_and_2_derivs[2], a_dt2quat,
                                   custom_approx);
    }
  }

  {
    INFO("QuaternionFunctionOfTime: Linear Omega");
    double t = 0.0;
    double expir_time = 0.5;
    // phi(t) = fac1 * t^2 + fac2 * t + fac3
    const double fac1 = 0.25;
    const double fac2 = 0.5;
    const double fac3 = 0.0;
    // omega(t) = 2*fac1 * t + fac2
    // dtomega(t) = 2 * fac1 (constant)
    // dt2omega(t) = 0.0

    DataVector init_quat{1.0, 0.0, 0.0, 0.0};
    DataVector init_angle{3, 0.0};
    DataVector init_omega{{0.0, 0.0, fac2}};
    DataVector init_dtomega{{0.0, 0.0, 2.0 * fac1}};
    // Construct QuaternionFunctionOfTime
    domain::FunctionsOfTime::QuaternionFunctionOfTime<2> qfot{
        t, std::array<DataVector, 1>{init_quat},
        std::array<DataVector, 3>{init_angle, init_omega, init_dtomega},
        expir_time};

    // Update internal PiecewisePolynomial with constant 2nd derivative so
    // omega is linear. This will automatically update the stored quaternions as
    // well
    const double time_step = 0.5;
    for (int i = 0; i < 15; i++) {
      t += time_step;
      expir_time += time_step;
      qfot.update(t, DataVector{{0.0, 0.0, 2.0 * fac1}}, expir_time);
      CHECK(qfot.expiration_after(t) == expir_time);
      CHECK(qfot.expiration_after(t + 0.5 * time_step) == expir_time);
      CHECK(qfot.expiration_after(t - 0.5 * time_step) == t);
    }
    CHECK(qfot.expiration_after(0.0) == time_step);

    // Get the quaternion and 2 derivatives at a certain time.
    double check_time = 5.398;
    const std::array<DataVector, 3> quat_func_and_2_derivs =
        qfot.quat_func_and_2_derivs(check_time);
    const std::array<DataVector, 3> quat_func_and_2_derivs2 =
        qfot.func_and_2_derivs(check_time);
    for (size_t i = 0; i < 3; i++) {
      CHECK_ITERABLE_APPROX(gsl::at(quat_func_and_2_derivs, i),
                            gsl::at(quat_func_and_2_derivs2, i));
    }

    // phi(t) = fac1 * t^2 + fac2 * t + fac3
    const double phi =
        fac1 * check_time * check_time + fac2 * check_time + fac3;
    // omega(t) = 2*fac1 * t + fac2
    const double omega = 2 * fac1 * check_time + fac2;
    const double dtomega = 2 * fac1;

    DataVector a_quat{{cos(0.5 * phi), 0.0, 0.0, sin(0.5 * phi)}};
    DataVector a_dtquat{
        {-0.5 * a_quat[3] * omega, 0.0, 0.0, 0.5 * a_quat[0] * omega}};
    DataVector a_dt2quat{{-0.5 * (a_dtquat[3] * omega + a_quat[3] * dtomega),
                          0.0, 0.0,
                          0.5 * (a_dtquat[0] * omega + a_quat[0] * dtomega)}};

    // Compare analytic solution to numerical
    Approx custom_approx = Approx::custom().epsilon(5.0e-12).scale(1.0);
    {
      INFO("  Compare quaternion");
      CHECK_ITERABLE_CUSTOM_APPROX(quat_func_and_2_derivs[0], a_quat,
                                   custom_approx);
    }
    {
      INFO("  Compare derivative of quaternion");
      CHECK_ITERABLE_CUSTOM_APPROX(quat_func_and_2_derivs[1], a_dtquat,
                                   custom_approx);
    }
    {
      INFO("  Compare second derivative of quaternion");
      CHECK_ITERABLE_CUSTOM_APPROX(quat_func_and_2_derivs[2], a_dt2quat,
                                   custom_approx);
    }
  }

  {
    INFO("QuaternionFunctionOfTime: Quadratic Omega");
    double t = 0.0;
    double expir_time = 0.5;
    // phi(t) = fac1 * t^3 + fac2 * t^2 + fac3 * t + fac4;
    const double fac1 = 0.2;
    const double fac2 = 0.3;
    const double fac3 = 0.4;
    const double fac4 = 0.0;
    // omega(t) = 3*fac1 * t^2 + 2*fac2 * t + fac3
    // dtomega(t) = 6*fac1 * t + 2* fac2
    // dt2omega(t) = 6 * fac1 (constant)

    DataVector init_quat{1.0, 0.0, 0.0, 0.0};
    DataVector init_angle{3, 0.0};
    DataVector init_omega{{0.0, 0.0, fac3}};
    DataVector init_dtomega{{0.0, 0.0, 2.0 * fac2}};
    DataVector init_dt2omega{{0.0, 0.0, 6.0 * fac1}};
    // Construct QuaternionFunctionOfTime
    domain::FunctionsOfTime::QuaternionFunctionOfTime<3> qfot{
        t, std::array<DataVector, 1>{init_quat},
        std::array<DataVector, 4>{init_angle, init_omega, init_dtomega,
                                  init_dt2omega},
        expir_time};

    // Update PiecewisePolynomial with constant 3rd derivative so omega is
    // quadratic. This will automatically update the stored quaternions as
    // well
    const double time_step = 0.5;
    for (int i = 0; i < 15; i++) {
      t += time_step;
      expir_time += time_step;
      qfot.update(t, DataVector{{0.0, 0.0, 6.0 * fac1}}, expir_time);
      CHECK(qfot.expiration_after(t) == expir_time);
      CHECK(qfot.expiration_after(t + 0.5 * time_step) == expir_time);
      CHECK(qfot.expiration_after(t - 0.5 * time_step) == t);
    }
    CHECK(qfot.expiration_after(0.0) == time_step);

    // Get the quaternion and 3 derivatives at a certain time.
    double check_time = 5.398;
    const std::array<DataVector, 3> quat_func_and_2_derivs =
        qfot.quat_func_and_2_derivs(check_time);
    const std::array<DataVector, 3> quat_func_and_2_derivs2 =
        qfot.func_and_2_derivs(check_time);
    for (size_t i = 0; i < 3; i++) {
      CHECK_ITERABLE_APPROX(gsl::at(quat_func_and_2_derivs, i),
                            gsl::at(quat_func_and_2_derivs2, i));
    }

    // phi(t) = fac1 * t^3 + fac2 * t^2 + fac3 * t + fac4;
    const double phi = fac1 * check_time * check_time * check_time +
                       fac2 * check_time * check_time + fac3 * check_time +
                       fac4;
    // omega(t) = 3*fac1 * t^2 + 2*fac2 * t + fac3
    // dtomega(t) = 6*fac1 * t + 2* fac2
    // dt2omega(t) = 6 * fac1 (constant)
    const double omega =
        3 * fac1 * check_time * check_time + 2 * fac2 * check_time + fac3;
    const double dtomega = 6 * fac1 * check_time + 2 * fac2;
    const double dt2omega = 6.0 * fac1;

    DataVector a_quat{{cos(0.5 * phi), 0.0, 0.0, sin(0.5 * phi)}};
    DataVector a_dtquat{
        {-0.5 * a_quat[3] * omega, 0.0, 0.0, 0.5 * a_quat[0] * omega}};
    DataVector a_dt2quat{{-0.5 * (a_dtquat[3] * omega + a_quat[3] * dtomega),
                          0.0, 0.0,
                          0.5 * (a_dtquat[0] * omega + a_quat[0] * dtomega)}};
    DataVector a_dt3quat{
        {-0.5 * (a_dt2quat[3] * omega + 2.0 * a_dtquat[3] * dtomega +
                 a_quat[3] * dt2omega),
         0.0, 0.0,
         0.5 * (a_dt2quat[0] * omega + 2.0 * a_dtquat[0] * dtomega +
                a_quat[0] * dt2omega)}};

    // Compare analytic solution to numerical
    Approx custom_approx = Approx::custom().epsilon(5.0e-12).scale(1.0);
    {
      INFO("  Compare quaternion");
      CHECK_ITERABLE_CUSTOM_APPROX(quat_func_and_2_derivs[0], a_quat,
                                   custom_approx);
    }
    {
      INFO("  Compare derivative of quaternion");
      CHECK_ITERABLE_CUSTOM_APPROX(quat_func_and_2_derivs[1], a_dtquat,
                                   custom_approx);
    }
    {
      INFO("  Compare second derivative of quaternion");
      CHECK_ITERABLE_CUSTOM_APPROX(quat_func_and_2_derivs[2], a_dt2quat,
                                   custom_approx);
    }
    {
      INFO("  Compare third derivative of quaternion");
      const std::array<DataVector, 4> quat_func_and_3_derivs =
          qfot.quat_func_and_3_derivs(check_time);
      CHECK_ITERABLE_CUSTOM_APPROX(quat_func_and_3_derivs[3], a_dt3quat,
                                   custom_approx);
    }
  }
  {
    INFO("QuaternionFunctionOfTime: No updates");
    const double initial_time = 0.0;
    const double final_time = 100.0;

    DataVector init_quat{1.0, 0.0, 0.0, 0.0};
    // We use zero because we are only concerned about checking the quaternion
    // at late times when it hasn't been updated
    DataVector three_zero{3, 0.0};
    // Construct QuaternionFunctionOfTime
    domain::FunctionsOfTime::QuaternionFunctionOfTime<3> qfot{
        initial_time, std::array<DataVector, 1>{init_quat},
        std::array<DataVector, 4>{three_zero, three_zero, three_zero,
                                  three_zero},
        std::numeric_limits<double>::infinity()};

    const std::array<DataVector, 3> quat_func_and_2_derivs =
        qfot.func_and_2_derivs(final_time);

    DataVector four_zero{4, 0.0};
    {
      INFO("  Compare quaternion");
      CHECK_ITERABLE_APPROX(quat_func_and_2_derivs[0], init_quat);
    }
    {
      INFO("  Compare derivative of quaternion");
      CHECK_ITERABLE_APPROX(quat_func_and_2_derivs[1], four_zero);
    }
    {
      INFO("  Compare second derivative of quaternion");
      CHECK_ITERABLE_APPROX(quat_func_and_2_derivs[2], four_zero);
    }
  }

  test_serialization_versioning();
  test_out_of_order_update();
}
