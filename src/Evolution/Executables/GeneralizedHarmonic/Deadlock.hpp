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
#include "Utilities/FileSystem.hpp"
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
  const std::string deadlock_dir{"deadlock"};
  if (file_system::check_if_dir_exists(deadlock_dir)) {
    file_system::rm(deadlock_dir, true);
  }

  file_system::create_directory(deadlock_dir);

  Parallel::simple_action<::deadlock::PrintFunctionsOfTime>(
      Parallel::get_parallel_component<
          observers::ObserverWriter<Metavariables>>(cache),
      deadlock_dir + "functions_of_time.out");

  if (alg::count(deadlocked_components, pretty_type::name<DgElementArray>()) ==
      1) {
    tmpl::for_each<ControlComponents>([&cache,
                                       &deadlock_dir](auto component_v) {
      using component = tmpl::type_from<decltype(component_v)>;
      Parallel::simple_action<control_system::Actions::PrintCurrentMeasurement>(
          Parallel::get_parallel_component<component>(cache),
          deadlock_dir + "/control_systems.out");
    });

    const std::string element_array_file =
        deadlock_dir + "/dg_element_array.out";
    if constexpr (Parallel::is_dg_element_collection_v<DgElementArray>) {
      Parallel::threaded_action<Parallel::Actions::SimpleActionOnElement<
          ::deadlock::PrintElementInfo, true>>(
          Parallel::get_parallel_component<DgElementArray>(cache),
          element_array_file);
    } else {
      Parallel::simple_action<::deadlock::PrintElementInfo>(
          Parallel::get_parallel_component<DgElementArray>(cache),
          element_array_file);
    }
  }
}
}  // namespace gh::deadlock
