// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <array>
#include <limits>
#include <memory>
#include <string>
#include <utility>

#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "Domain/Creators/Tags/FunctionsOfTime.hpp"
#include "Domain/FunctionsOfTime/PiecewisePolynomial.hpp"
#include "Domain/FunctionsOfTime/Tags.hpp"
#include "Framework/ActionTesting.hpp"
#include "Framework/MockDistributedObject.hpp"
#include "Framework/MockRuntimeSystem.hpp"
#include "Framework/MockRuntimeSystemFreeFunctions.hpp"
#include "Parallel/Callback.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "ParallelAlgorithms/Actions/FunctionsOfTimeAreReady.hpp"
#include "Time/Tags/Time.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/TMPL.hpp"

namespace {
using FunctionMap = domain::Tags::FunctionsOfTimeInitialize::type;

struct OtherFunctionsOfTime : db::SimpleTag {
  using type = FunctionMap;
};

struct UpdateFoT {
  static void apply(const gsl::not_null<FunctionMap*> functions,
                    const std::string& name, const double expiration) {
    const double current_expiration = functions->at(name)->time_bounds()[1];
    // Update value doesn't matter
    (*functions)
        .at(name)
        ->update(current_expiration, DataVector{0.0}, expiration);
  }
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
size_t simple_action_no_args_call_count = 0;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
size_t threaded_action_no_args_call_count = 0;

struct SimpleActionNoArgs {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex>
  static void apply(db::DataBox<DbTags>& /*box*/,
                    Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/) {
    ++simple_action_no_args_call_count;
  }
};

struct ThreadedActionNoArgs {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex>
  static void apply(db::DataBox<DbTags>& /*box*/,
                    Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/,
                    const gsl::not_null<Parallel::NodeLock*> /*node_lock*/) {
    ++threaded_action_no_args_call_count;
  }
};

// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
size_t simple_action_args_call_count = 0;
// NOLINTNEXTLINE(cppcoreguidelines-avoid-non-const-global-variables)
size_t threaded_action_args_call_count = 0;

struct SimpleActionArgs {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex>
  static void apply(db::DataBox<DbTags>& /*box*/,
                    Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/, const size_t size,
                    const DataVector& some_data) {
    ++simple_action_args_call_count;
    CHECK(some_data.size() == size);
  }
};

struct ThreadedActionArgs {
  template <typename ParallelComponent, typename DbTags, typename Metavariables,
            typename ArrayIndex>
  static void apply(db::DataBox<DbTags>& /*box*/,
                    Parallel::GlobalCache<Metavariables>& /*cache*/,
                    const ArrayIndex& /*array_index*/,
                    const gsl::not_null<Parallel::NodeLock*> /*node_lock*/,
                    const size_t size, const DataVector& some_data) {
    ++threaded_action_args_call_count;
    CHECK(some_data.size() == size);
  }
};

template <typename Metavariables, size_t Index>
struct Component {
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = int;

  using simple_tags_from_options = tmpl::list<Tags::Time>;
  using mutable_global_cache_tags =
      tmpl::list<domain::Tags::FunctionsOfTimeInitialize, OtherFunctionsOfTime>;

  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      Parallel::Phase::Testing,
      tmpl::list<domain::Actions::CheckFunctionsOfTimeAreReady<3>>>>;
};

template <typename Metavariables, size_t Index>
struct EmptyComponent : Component<Metavariables, Index> {
  using mutable_global_cache_tags = tmpl::list<>;
};
struct EmptyMetavars {
  using component_list = tmpl::list<EmptyComponent<EmptyMetavars, 0_st>>;
};

struct Metavariables {
  using component_list = tmpl::list<Component<Metavariables, 0_st>,
                                    Component<Metavariables, 1_st>>;
};
}  // namespace

SPECTRE_TEST_CASE("Unit.Domain.FunctionsOfTimeAreReady", "[Domain][Unit]") {
  register_classes_with_charm<
      domain::FunctionsOfTime::PiecewisePolynomial<2>>();
  using MockRuntimeSystem = ActionTesting::MockRuntimeSystem<Metavariables>;
  using component0 = Component<Metavariables, 0_st>;
  using component1 = Component<Metavariables, 1_st>;
  const component0* const component_0_p = nullptr;
  const component1* const component_1_p = nullptr;
  // Not used when not using nodegroup element components
  constexpr size_t dim = std::numeric_limits<size_t>::max();

  const std::array<DataVector, 3> fot_init{{{0.0}, {0.0}, {0.0}}};
  FunctionMap functions_of_time{};
  functions_of_time["FunctionA"] =
      std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<2>>(
          0.0, fot_init, 1.0);
  functions_of_time["FunctionB"] =
      std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<2>>(
          0.0, fot_init, 1.0);

  FunctionMap other_functions_of_time{};
  other_functions_of_time["OtherA"] =
      std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<2>>(
          0.0, fot_init, 0.1);
  other_functions_of_time["OtherB"] =
      std::make_unique<domain::FunctionsOfTime::PiecewisePolynomial<2>>(
          0.0, fot_init, 0.1);

  MockRuntimeSystem runner{
      {}, {std::move(functions_of_time), std::move(other_functions_of_time)}};
  ActionTesting::emplace_array_component<component0>(
      make_not_null(&runner), ActionTesting::NodeId{0},
      ActionTesting::LocalCoreId{0}, 0, 2.0);
  ActionTesting::emplace_array_component<component1>(
      make_not_null(&runner), ActionTesting::NodeId{0},
      ActionTesting::LocalCoreId{0}, 0, 2.0);
  ActionTesting::set_phase(make_not_null(&runner), Parallel::Phase::Testing);

  auto& cache = ActionTesting::cache<component0>(runner, 0);

  // Test the algorithm callback free function.  For this section of the test,
  // the "standard" tags domain::Tags::FunctionsOfTime is ready, but this code
  // should never examine it.
  {
    INFO("Testing perform algorithm callback support");
    // Neither function ready
    CHECK(not domain::functions_of_time_are_ready_algorithm_callback<
          OtherFunctionsOfTime, dim>(cache, 0, component_0_p, 0.5,
                                     std::nullopt));
    CHECK(not domain::functions_of_time_are_ready_algorithm_callback<
          OtherFunctionsOfTime, dim>(cache, 0, component_0_p, 0.5,
                                     std::unordered_set{"OtherA"s}));

    // Make OtherA ready
    Parallel::mutate<OtherFunctionsOfTime, UpdateFoT>(cache, "OtherA"s, 123.0);

    CHECK(domain::functions_of_time_are_ready_algorithm_callback<
          OtherFunctionsOfTime, dim>(cache, 0, component_0_p, 0.5,
                                     std::unordered_set{"OtherA"s}));
    CHECK(not domain::functions_of_time_are_ready_algorithm_callback<
          OtherFunctionsOfTime, dim>(cache, 0, component_0_p, 0.5,
                                     std::unordered_set{"OtherA"s, "OtherB"s}));

    // Make OtherB ready
    Parallel::mutate<OtherFunctionsOfTime, UpdateFoT>(cache, "OtherB"s, 456.0);

    CHECK(domain::functions_of_time_are_ready_algorithm_callback<
          OtherFunctionsOfTime, dim>(cache, 0, component_0_p, 0.5,
                                     std::unordered_set{"OtherA"s, "OtherB"s}));
    CHECK(domain::functions_of_time_are_ready_algorithm_callback<
          OtherFunctionsOfTime, dim>(cache, 0, component_0_p, 0.5,
                                     std::nullopt));
  }

  // Test the action.  This should automatically look at
  // domain::Tags::FunctionsOfTime.
  {
    // Neither function ready
    CHECK(not ActionTesting::next_action_if_ready<component0>(
        make_not_null(&runner), 0));

    // Make OtherA ready
    Parallel::mutate<domain::Tags::FunctionsOfTime, UpdateFoT>(
        cache, "FunctionA"s, 5.0);

    CHECK(not ActionTesting::next_action_if_ready<component0>(
        make_not_null(&runner), 0));

    // Make OtherB ready
    Parallel::mutate<domain::Tags::FunctionsOfTime, UpdateFoT>(
        cache, "FunctionB"s, 10.0);

    CHECK(ActionTesting::next_action_if_ready<component0>(
        make_not_null(&runner), 0));
  }

  // Test simple action callback free function
  {
    INFO("Testing simple action support");
    DataVector data{6, 0.0};
    const size_t size = data.size();

    // No callbacks should be registered
    CHECK(domain::functions_of_time_are_ready_simple_action_callback<
          domain::Tags::FunctionsOfTime, SimpleActionNoArgs>(
        cache, 0, component_0_p, 5.0, std::nullopt));
    CHECK(ActionTesting::number_of_queued_simple_actions<component0>(runner,
                                                                     0) == 0);
    CHECK(simple_action_no_args_call_count == 0_st);
    CHECK(domain::functions_of_time_are_ready_simple_action_callback<
          domain::Tags::FunctionsOfTime, SimpleActionArgs>(
        cache, 0, component_0_p, 5.0, std::nullopt, size, data));
    CHECK(ActionTesting::number_of_queued_simple_actions<component0>(runner,
                                                                     0) == 0);
    CHECK(simple_action_args_call_count == 0_st);

    // Again no callbacks should be registered
    CHECK(domain::functions_of_time_are_ready_simple_action_callback<
          domain::Tags::FunctionsOfTime, SimpleActionNoArgs>(
        cache, 0, component_0_p, 6.0, std::unordered_set{"FunctionB"s}));
    CHECK(ActionTesting::number_of_queued_simple_actions<component0>(runner,
                                                                     0) == 0);
    CHECK(simple_action_no_args_call_count == 0_st);
    CHECK(domain::functions_of_time_are_ready_simple_action_callback<
          domain::Tags::FunctionsOfTime, SimpleActionArgs>(
        cache, 0, component_0_p, 6.0, std::unordered_set{"FunctionB"s}, size,
        data));
    CHECK(ActionTesting::number_of_queued_simple_actions<component0>(runner,
                                                                     0) == 0);
    CHECK(simple_action_args_call_count == 0_st);

    // Evaluate at time when A isn't ready. Can't have two different
    // callbacks on same component so we use different components
    CHECK(not domain::functions_of_time_are_ready_simple_action_callback<
          domain::Tags::FunctionsOfTime, SimpleActionNoArgs>(
        cache, 0, component_0_p, 6.0, std::unordered_set{"FunctionA"s}));
    CHECK(ActionTesting::number_of_queued_simple_actions<component0>(runner,
                                                                     0) == 0);
    CHECK(simple_action_no_args_call_count == 0_st);
    CHECK(not domain::functions_of_time_are_ready_simple_action_callback<
          domain::Tags::FunctionsOfTime, SimpleActionArgs>(
        cache, 0, component_1_p, 6.0, std::unordered_set{"FunctionA"s}, size,
        data));
    CHECK(ActionTesting::number_of_queued_simple_actions<component1>(runner,
                                                                     0) == 0);
    CHECK(simple_action_args_call_count == 0_st);

    // Make FunctionA valid again
    Parallel::mutate<domain::Tags::FunctionsOfTime, UpdateFoT>(
        cache, "FunctionA"s, 10.0);
    // Both actions should have been queued
    CHECK(ActionTesting::number_of_queued_simple_actions<component0>(runner,
                                                                     0) == 1);
    CHECK(ActionTesting::number_of_queued_simple_actions<component1>(runner,
                                                                     0) == 1);
    CHECK(simple_action_no_args_call_count == 0_st);
    CHECK(simple_action_args_call_count == 0_st);
    ActionTesting::invoke_queued_simple_action<component0>(
        make_not_null(&runner), 0);
    // Only one ran
    CHECK(simple_action_no_args_call_count == 1_st);
    CHECK(simple_action_args_call_count == 0_st);
    ActionTesting::invoke_queued_simple_action<component1>(
        make_not_null(&runner), 0);
    // Both ran
    CHECK(simple_action_no_args_call_count == 1_st);
    CHECK(simple_action_args_call_count == 1_st);
    CHECK(domain::functions_of_time_are_ready_simple_action_callback<
          domain::Tags::FunctionsOfTime, SimpleActionNoArgs>(
        cache, 0, component_0_p, 6.0, std::unordered_set{"FunctionA"s}));
    CHECK(ActionTesting::number_of_queued_simple_actions<component0>(runner,
                                                                     0) == 0);

    // Evaluate at time when neither are ready.
    CHECK(not domain::functions_of_time_are_ready_simple_action_callback<
          domain::Tags::FunctionsOfTime, SimpleActionNoArgs>(
        cache, 0, component_0_p, 11.0, std::nullopt));
    CHECK(ActionTesting::number_of_queued_simple_actions<component0>(runner,
                                                                     0) == 0);
    CHECK(simple_action_no_args_call_count == 1_st);
    CHECK(not domain::functions_of_time_are_ready_simple_action_callback<
          domain::Tags::FunctionsOfTime, SimpleActionArgs>(
        cache, 0, component_1_p, 11.0, std::nullopt, size, data));
    CHECK(ActionTesting::number_of_queued_simple_actions<component1>(runner,
                                                                     0) == 0);
    CHECK(simple_action_args_call_count == 1_st);

    // Make A valid
    Parallel::mutate<domain::Tags::FunctionsOfTime, UpdateFoT>(
        cache, "FunctionA"s, 15.0);
    // Both actions should have been queued
    CHECK(ActionTesting::number_of_queued_simple_actions<component0>(runner,
                                                                     0) == 1);
    CHECK(ActionTesting::number_of_queued_simple_actions<component1>(runner,
                                                                     0) == 1);
    CHECK(simple_action_no_args_call_count == 1_st);
    CHECK(simple_action_args_call_count == 1_st);
    ActionTesting::invoke_queued_simple_action<component0>(
        make_not_null(&runner), 0);
    CHECK(simple_action_no_args_call_count == 2_st);
    CHECK(simple_action_args_call_count == 1_st);
    ActionTesting::invoke_queued_simple_action<component1>(
        make_not_null(&runner), 0);
    CHECK(simple_action_no_args_call_count == 2_st);
    CHECK(simple_action_args_call_count == 2_st);

    // Make B valid. Nothing should have happened
    Parallel::mutate<domain::Tags::FunctionsOfTime, UpdateFoT>(
        cache, "FunctionB"s, 15.0);
    CHECK(ActionTesting::number_of_queued_simple_actions<component0>(runner,
                                                                     0) == 0);
    CHECK(ActionTesting::number_of_queued_simple_actions<component1>(runner,
                                                                     0) == 0);
    // No actions should have run
    CHECK(simple_action_no_args_call_count == 2_st);
    CHECK(simple_action_args_call_count == 2_st);
  }

  // Test threaded action callback free function
  {
    INFO("Testing threaded action support");
    DataVector data{6, 0.0};
    const size_t size = data.size();

    // No callbacks should be registered
    CHECK(domain::functions_of_time_are_ready_threaded_action_callback<
          domain::Tags::FunctionsOfTime, ThreadedActionNoArgs>(
        cache, 0, component_0_p, 5.0, std::nullopt));
    CHECK(ActionTesting::number_of_queued_threaded_actions<component0>(runner,
                                                                       0) == 0);
    CHECK(threaded_action_no_args_call_count == 0_st);
    CHECK(domain::functions_of_time_are_ready_threaded_action_callback<
          domain::Tags::FunctionsOfTime, ThreadedActionArgs>(
        cache, 0, component_0_p, 5.0, std::nullopt, size, data));
    CHECK(ActionTesting::number_of_queued_threaded_actions<component0>(runner,
                                                                       0) == 0);
    CHECK(threaded_action_args_call_count == 0_st);

    // Again no callbacks should be registered
    CHECK(domain::functions_of_time_are_ready_threaded_action_callback<
          domain::Tags::FunctionsOfTime, ThreadedActionNoArgs>(
        cache, 0, component_0_p, 6.0, std::unordered_set{"FunctionB"s}));
    CHECK(ActionTesting::number_of_queued_threaded_actions<component0>(runner,
                                                                       0) == 0);
    CHECK(threaded_action_no_args_call_count == 0_st);
    CHECK(domain::functions_of_time_are_ready_threaded_action_callback<
          domain::Tags::FunctionsOfTime, ThreadedActionArgs>(
        cache, 0, component_0_p, 6.0, std::unordered_set{"FunctionB"s}, size,
        data));
    CHECK(ActionTesting::number_of_queued_threaded_actions<component0>(runner,
                                                                       0) == 0);
    CHECK(threaded_action_args_call_count == 0_st);

    // Evaluate at time when A isn't ready. Can't have two different
    // callbacks on same component so we use different components
    CHECK(not domain::functions_of_time_are_ready_threaded_action_callback<
          domain::Tags::FunctionsOfTime, ThreadedActionNoArgs>(
        cache, 0, component_0_p, 16.0, std::unordered_set{"FunctionA"s}));
    CHECK(ActionTesting::number_of_queued_threaded_actions<component0>(runner,
                                                                       0) == 0);
    CHECK(threaded_action_no_args_call_count == 0_st);
    CHECK(not domain::functions_of_time_are_ready_threaded_action_callback<
          domain::Tags::FunctionsOfTime, ThreadedActionArgs>(
        cache, 0, component_1_p, 16.0, std::unordered_set{"FunctionA"s}, size,
        data));
    CHECK(ActionTesting::number_of_queued_threaded_actions<component1>(runner,
                                                                       0) == 0);
    CHECK(threaded_action_args_call_count == 0_st);

    // Make FunctionA valid again
    Parallel::mutate<domain::Tags::FunctionsOfTime, UpdateFoT>(
        cache, "FunctionA"s, 20.0);
    // Both actions should have been queued
    CHECK(ActionTesting::number_of_queued_threaded_actions<component0>(runner,
                                                                       0) == 1);
    CHECK(ActionTesting::number_of_queued_threaded_actions<component1>(runner,
                                                                       0) == 1);
    CHECK(threaded_action_no_args_call_count == 0_st);
    CHECK(threaded_action_args_call_count == 0_st);
    ActionTesting::invoke_queued_threaded_action<component0>(
        make_not_null(&runner), 0);
    // Only one ran
    CHECK(threaded_action_no_args_call_count == 1_st);
    CHECK(threaded_action_args_call_count == 0_st);
    ActionTesting::invoke_queued_threaded_action<component1>(
        make_not_null(&runner), 0);
    // Both ran
    CHECK(threaded_action_no_args_call_count == 1_st);
    CHECK(threaded_action_args_call_count == 1_st);
    CHECK(domain::functions_of_time_are_ready_threaded_action_callback<
          domain::Tags::FunctionsOfTime, ThreadedActionNoArgs>(
        cache, 0, component_0_p, 16.0, std::unordered_set{"FunctionA"s}));
    CHECK(ActionTesting::number_of_queued_threaded_actions<component0>(runner,
                                                                       0) == 0);

    // Evaluate at time when neither are ready.
    CHECK(not domain::functions_of_time_are_ready_threaded_action_callback<
          domain::Tags::FunctionsOfTime, ThreadedActionNoArgs>(
        cache, 0, component_0_p, 21.0, std::nullopt));
    CHECK(ActionTesting::number_of_queued_threaded_actions<component0>(runner,
                                                                       0) == 0);
    CHECK(threaded_action_no_args_call_count == 1_st);
    CHECK(not domain::functions_of_time_are_ready_threaded_action_callback<
          domain::Tags::FunctionsOfTime, ThreadedActionArgs>(
        cache, 0, component_1_p, 21.0, std::nullopt, size, data));
    CHECK(ActionTesting::number_of_queued_threaded_actions<component1>(runner,
                                                                       0) == 0);
    CHECK(threaded_action_args_call_count == 1_st);

    // Make A valid
    Parallel::mutate<domain::Tags::FunctionsOfTime, UpdateFoT>(
        cache, "FunctionA"s, 25.0);
    // Both actions should have been queued
    CHECK(ActionTesting::number_of_queued_threaded_actions<component0>(runner,
                                                                       0) == 1);
    CHECK(ActionTesting::number_of_queued_threaded_actions<component1>(runner,
                                                                       0) == 1);
    CHECK(threaded_action_no_args_call_count == 1_st);
    CHECK(threaded_action_args_call_count == 1_st);
    ActionTesting::invoke_queued_threaded_action<component0>(
        make_not_null(&runner), 0);
    CHECK(threaded_action_no_args_call_count == 2_st);
    CHECK(threaded_action_args_call_count == 1_st);
    ActionTesting::invoke_queued_threaded_action<component1>(
        make_not_null(&runner), 0);
    CHECK(threaded_action_no_args_call_count == 2_st);
    CHECK(threaded_action_args_call_count == 2_st);

    // Make B valid. Nothing should have happened
    Parallel::mutate<domain::Tags::FunctionsOfTime, UpdateFoT>(
        cache, "FunctionB"s, 25.0);
    CHECK(ActionTesting::number_of_queued_threaded_actions<component0>(runner,
                                                                       0) == 0);
    CHECK(ActionTesting::number_of_queued_threaded_actions<component1>(runner,
                                                                       0) == 0);
    // No actions should have run
    CHECK(threaded_action_no_args_call_count == 2_st);
    CHECK(threaded_action_args_call_count == 2_st);
  }

  // No FoTs in the cache
  {
    using EmptyMockRuntimeSystem =
        ActionTesting::MockRuntimeSystem<EmptyMetavars>;
    using empty_comp = EmptyComponent<EmptyMetavars, 0_st>;
    const empty_comp* const empty_component_p = nullptr;

    EmptyMockRuntimeSystem empty_runner{{}};
    ActionTesting::emplace_array_component<empty_comp>(
        make_not_null(&empty_runner), ActionTesting::NodeId{0},
        ActionTesting::LocalCoreId{0}, 0, 2.0);

    auto& empty_cache = ActionTesting::cache<empty_comp>(empty_runner, 0);
    CHECK(domain::functions_of_time_are_ready_algorithm_callback<
          domain::Tags::FunctionsOfTime, dim>(empty_cache, 0, empty_component_p,
                                              6.0));
    CHECK(domain::functions_of_time_are_ready_simple_action_callback<
          domain::Tags::FunctionsOfTime, SimpleActionNoArgs>(
        empty_cache, 0, empty_component_p, 6.0, std::nullopt));
    CHECK(domain::functions_of_time_are_ready_simple_action_callback<
          domain::Tags::FunctionsOfTime, SimpleActionArgs>(
        empty_cache, 0, empty_component_p, 6.0, std::nullopt, 0_st,
        DataVector{}));
  }
}
