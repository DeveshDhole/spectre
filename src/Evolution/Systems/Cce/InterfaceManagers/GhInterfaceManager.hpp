// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <memory>
#include <optional>
#include <tuple>

#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Tensor/TypeAliases.hpp"
#include "Evolution/Systems/GeneralizedHarmonic/Tags.hpp"
#include "Options/String.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Time/TimeStepId.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"

namespace Cce {
/// \brief Code for interfacing between the characteristic and Cachy systems.
namespace InterfaceManagers {

/// \cond
class GhLocalTimeStepping;
class GhLockstep;
/// \endcond

/*!
 * \brief Abstract base class for storage and retrieval of generalized harmonic
 * quantities communicated from a Cauchy simulation to the Cce system.
 *
 * \details The functions that are required to be overriden in the derived
 * classes are:
 * - `GhInterfaceManager::get_clone()`: should return a
 * `std::unique_ptr<GhInterfaceManager>` with cloned state.
 * - `GhInterfaceManager::insert_gh_data()`: should store the portions
 * of the provided generalized harmonic data that are required to provide useful
 * boundary values for the CCE evolution at requested timesteps.
 * - `GhInterfaceManager::request_gh_data()`: should register requests
 * from the CCE evolution for boundary data.
 * - `GhInterfaceManager::retrieve_and_remove_first_ready_gh_data()`:
 * should return a `std::optional<std::tuple<TimeStepId,
 * tnsr::aa<DataVector, 3>, tnsr::iaa<DataVector, 3>, tnsr::aa<DataVector, 3>>>`
 * containing the boundary data associated with the oldest requested timestep if
 * enough data has been supplied via `insert_gh_data()` to determine the
 * boundary data. Otherwise, return a `std::nullopt` to indicate that the CCE
 * system must continue waiting for generalized harmonic input.
 * - `GhInterfaceManager::number_of_pending_requests()`: should return
 * the number of requests that have been registered to the class that do not yet
 * been retrieved via `retrieve_and_remove_first_ready_gh_data()`.
 * - `GhInterfaceManager::number_of_gh_times()`: should return the
 * number of time steps sent to `insert_gh_data()` that have not yet been
 * retrieved via `retrieve_and_remove_first_ready_gh_data()`.
 */
class GhInterfaceManager : public PUP::able {
 public:
  using gh_variables = Variables<
      tmpl::list<gr::Tags::SpacetimeMetric<DataVector, 3>,
                 gh::Tags::Pi<DataVector, 3>, gh::Tags::Phi<DataVector, 3>>>;

  WRAPPED_PUPable_abstract(GhInterfaceManager);  // NOLINT

  virtual std::unique_ptr<GhInterfaceManager> get_clone() const = 0;

  virtual void request_gh_data(const TimeStepId&) = 0;

  virtual auto retrieve_and_remove_first_ready_gh_data()
      -> std::optional<std::tuple<TimeStepId, gh_variables>> = 0;

  virtual size_t number_of_pending_requests() const = 0;

  virtual size_t number_of_gh_times() const = 0;
};

}  // namespace InterfaceManagers
}  // namespace Cce
