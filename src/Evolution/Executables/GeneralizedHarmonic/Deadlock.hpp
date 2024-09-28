// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <string>
#include <vector>

#include "ControlSystem/Actions/PrintCurrentMeasurement.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Evolution/Deadlock/PrintDgElementArray.hpp"
#include "Evolution/Deadlock/PrintFunctionsOfTime.hpp"
#include "Parallel/ArrayCollection/IsDgElementCollection.hpp"
#include "Parallel/ArrayCollection/SimpleActionOnElement.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Invoke.hpp"
#include "Utilities/PrettyType.hpp"
#include "Utilities/TMPL.hpp"

/// \cond
namespace observers {
template <class Metavariables>
struct ObserverWriter;
}  // namespace observers
/// \endcond

namespace gh::deadlock {
template <typename DgElementArray, typename ControlComponents,
          typename Metavariables>
void run_deadlock_analysis_simple_actions(
    Parallel::GlobalCache<Metavariables>& cache,
    const std::vector<std::string>& deadlocked_components) {
  Parallel::simple_action<::deadlock::PrintFunctionsOfTime>(
      Parallel::get_parallel_component<
          observers::ObserverWriter<Metavariables>>(cache));

  if (alg::count(deadlocked_components, pretty_type::name<DgElementArray>()) ==
      1) {
    tmpl::for_each<ControlComponents>([&cache](auto component_v) {
      using component = tmpl::type_from<decltype(component_v)>;
      Parallel::simple_action<control_system::Actions::PrintCurrentMeasurement>(
          Parallel::get_parallel_component<component>(cache));
    });

    if constexpr (Parallel::is_dg_element_collection_v<DgElementArray>) {
      Parallel::threaded_action<Parallel::Actions::SimpleActionOnElement<
          ::deadlock::PrintElementInfo, true>>(
          Parallel::get_parallel_component<DgElementArray>(cache));
    } else {
      Parallel::simple_action<::deadlock::PrintElementInfo>(
          Parallel::get_parallel_component<DgElementArray>(cache));
    }
  }
}
}  // namespace gh::deadlock
