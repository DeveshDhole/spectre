// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Helpers/PointwiseFunctions/PostNewtonian/BinaryTrajectories.hpp"

#include <array>
#include <cmath>
#include <utility>

#include "Utilities/ConstantExpressions.hpp"

BinaryTrajectories::BinaryTrajectories(
    double initial_separation,
    const std::array<double, 3>& center_of_mass_velocity, bool newtonian)
    : initial_separation_fourth_power_{square(square(initial_separation))},
      center_of_mass_velocity_(center_of_mass_velocity),
      newtonian_(newtonian) {}

double BinaryTrajectories::separation(const double time) const {
  const double pn_correction_term = newtonian_ ? 0.0 : 12.8 * time;
  return pow(initial_separation_fourth_power_ - pn_correction_term, 0.25);
}

double BinaryTrajectories::orbital_frequency(const double time) const {
  return pow(separation(time), -1.5);
}

double BinaryTrajectories::angular_velocity(const double time) const {
  // This is d/dt(orbital_frequency) if we are using PN, but 0 if it's newtonian
  const double pn_correction_term =
      newtonian_
          ? 0.0
          : 4.8 * pow(initial_separation_fourth_power_ - 12.8 * time, -1.375);
  return orbital_frequency(time) + pn_correction_term * time;
}

std::pair<std::array<double, 3>, std::array<double, 3>>
BinaryTrajectories::positions(const double time) const {
  const double sep = separation(time);
  return position_impl(time, sep);
}

std::pair<std::array<double, 3>, std::array<double, 3>>
BinaryTrajectories::positions_no_expansion(const double time) const {
  // Separation stays constant while orbital frequency follows PN (or newtonian)
  // values
  const double sep = pow(initial_separation_fourth_power_, 0.25);
  return position_impl(time, sep);
}

std::pair<std::array<double, 3>, std::array<double, 3>>
BinaryTrajectories::position_impl(const double time,
                                  const double separation) const {
  const double orbital_freq = orbital_frequency(time);
  double xA = 0.5 * separation * cos(orbital_freq * time);
  double yA = 0.5 * separation * sin(orbital_freq * time);
  double xB = -xA;
  double yB = -yA;
  xB += center_of_mass_velocity_[0] * time;
  xA += center_of_mass_velocity_[0] * time;
  yB += center_of_mass_velocity_[1] * time;
  yA += center_of_mass_velocity_[1] * time;
  const double zB = center_of_mass_velocity_[2] * time;
  const double zA = zB;
  return {{xA, yA, zA}, {xB, yB, zB}};
}
