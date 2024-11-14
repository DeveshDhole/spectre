// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <optional>
#include <pup.h>

#include "DataStructures/DataVector.hpp"
#include "Domain/CoordinateMaps/TimeDependent/ShapeMapTransitionFunctions/ShapeMapTransitionFunction.hpp"

namespace domain::CoordinateMaps::ShapeMapTransitionFunctions {

/*!
 * \ingroup CoordMapsTimeDependentGroup
 * \brief A transition function that falls off as $f(r) = ar + b$.
 *
 * \details The coefficients $a$ and $b$ are chosen so that the function $f(r) =
 * ar + b$ falls off linearly from 1 at `r_min` to 0 at `r_max`. The
 * coefficients are
 *
 * \f{align}{
 * a &= \frac{-1}{r_{\text{max}} - r_{\text{min}}} \\
 * b &= \frac{r_{\text{max}}}{r_{\text{max}} - r_{\text{min}}} = -a
 * r_{\text{max}}
 * \f}
 *
 * If \p reverse is set to `true`, then the function falls off from 0 at
 * `r_min` to 1 at `r_max`. To do this, the coefficients are modified as
 * $a \rightarrow -a$ and $b \rightarrow 1-b$.
 *
 * The function can be called beyond `r_min` and `r_max`. Within `r_min` the
 * value is 1, and outside `r_max` the value is 0. This is reversed if \p
 * reverse is true. However, the gradient function cannot be called with a point
 * beyond `r_min` and `r_max`.
 */
class SphereTransition final : public ShapeMapTransitionFunction {
 public:
  explicit SphereTransition() = default;
  SphereTransition(double r_min, double r_max, bool reverse = false);

  double operator()(const std::array<double, 3>& source_coords) const override;
  DataVector operator()(
      const std::array<DataVector, 3>& source_coords) const override;

  std::optional<double> original_radius_over_radius(
      const std::array<double, 3>& target_coords,
      double radial_distortion) const override;

  std::array<double, 3> gradient(
      const std::array<double, 3>& source_coords) const override;
  std::array<DataVector, 3> gradient(
      const std::array<DataVector, 3>& source_coords) const override;

  WRAPPED_PUPable_decl_template(SphereTransition);
  explicit SphereTransition(CkMigrateMessage* const msg);
  void pup(PUP::er& p) override;

  std::unique_ptr<ShapeMapTransitionFunction> get_clone() const override {
    return std::make_unique<SphereTransition>(*this);
  }

  bool operator==(const ShapeMapTransitionFunction& other) const override;
  bool operator!=(const ShapeMapTransitionFunction& other) const override;

 private:
  template <typename T>
  T call_impl(const std::array<T, 3>& source_coords) const;

  template <typename T>
  std::array<T, 3> gradient_impl(const std::array<T, 3>& source_coords) const;

  // checks that the magnitudes are all between `r_min_` and `r_max_`
  template <typename T>
  void check_magnitudes(const T& mag, bool check_bounds) const;

  double r_min_{};
  double r_max_{};
  double a_{};
  double b_{};
  static constexpr double eps_ = std::numeric_limits<double>::epsilon() * 100;
};
}  // namespace domain::CoordinateMaps::ShapeMapTransitionFunctions
