// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <cstddef>

#include "DataStructures/DataBox/DataBox.hpp"
#include "Domain/FunctionsOfTime/OutputTimeBounds.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Info.hpp"
#include "Parallel/Printf/Printf.hpp"

namespace deadlock {
/*!
 * \brief Simple action that will print the `domain::Tags::FunctionsOfTime` time
 * bounds for each node of a simulation.
 */
struct PrintFunctionsOfTime {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex>
  static void apply(db::DataBox<DbTags>& /*box*/,
                    Parallel::GlobalCache<Metavariables>& cache,
                    const ArrayIndex& /*array_index*/) {
    const auto& functions_of_time =
        Parallel::get<::domain::Tags::FunctionsOfTime>(cache);

    const std::string time_bounds =
        ::domain::FunctionsOfTime::output_time_bounds(functions_of_time);

    Parallel::printf("Node %zu\n%s\n", Parallel::my_node<size_t>(cache),
                     time_bounds);
  }
};
}  // namespace deadlock
