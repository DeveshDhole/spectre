// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <limits>
#include <memory>
#include <optional>
#include <pup.h>
#include <string>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "DataStructures/DataBox/DataBox.hpp"
#include "DataStructures/DataBox/DataBoxTag.hpp"
#include "DataStructures/DataBox/PrefixHelpers.hpp"
#include "DataStructures/DataBox/Prefixes.hpp"
#include "DataStructures/DataBox/Tag.hpp"
#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "DataStructures/Variables.hpp"
#include "DataStructures/VariablesTag.hpp"
#include "Domain/Structure/Direction.hpp"
#include "Domain/Structure/DirectionMap.hpp"
#include "Domain/Structure/DirectionalIdMap.hpp"
#include "Domain/Structure/Element.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Structure/Neighbors.hpp"
#include "Domain/Tags.hpp"
#include "Domain/Tags/NeighborMesh.hpp"
#include "Evolution/Actions/RunEventsAndDenseTriggers.hpp"
#include "Framework/ActionTesting.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Options/Protocols/FactoryCreation.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "ParallelAlgorithms/Amr/Protocols/Projector.hpp"
#include "ParallelAlgorithms/EventsAndDenseTriggers/DenseTrigger.hpp"
#include "ParallelAlgorithms/EventsAndDenseTriggers/EventsAndDenseTriggers.hpp"
#include "ParallelAlgorithms/EventsAndDenseTriggers/Tags.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "Time/History.hpp"
#include "Time/Slab.hpp"
#include "Time/Tags/HistoryEvolvedVariables.hpp"
#include "Time/Tags/Time.hpp"
#include "Time/Tags/TimeStep.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Time/Tags/TimeStepper.hpp"
#include "Time/Time.hpp"
#include "Time/TimeStepId.hpp"
#include "Time/TimeSteppers/AdamsBashforth.hpp"
#include "Time/TimeSteppers/Rk3HesthavenSsp.hpp"
#include "Time/TimeSteppers/TimeStepper.hpp"
#include "Utilities/GetOutput.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/MakeVector.hpp"
#include "Utilities/ProtocolHelpers.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"
#include "Utilities/Serialization/RegisterDerivedClassesWithCharm.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

namespace Parallel {
template <typename Metavariables>
class GlobalCache;
}  // namespace Parallel

namespace {
struct EvolvedVar : db::SimpleTag {
  using type = Scalar<DataVector>;
};

using EvolvedVariables = Variables<tmpl::list<EvolvedVar>>;

template <typename T, typename Label = void>
struct PostprocessedVar : db::SimpleTag {
  using type = T;
};

namespace labels {
struct A;
struct B;
}  // namespace labels

using extra_data =
    tmpl::list<PostprocessedVar<Scalar<DataVector>, labels::A>,
               PostprocessedVar<Scalar<DataVector>, labels::B>,
               PostprocessedVar<Scalar<double>>, PostprocessedVar<std::string>>;
using all_data = tmpl::push_front<extra_data, ::Tags::Time, EvolvedVar>;
using DataTuple = tuples::tagged_tuple_from_typelist<all_data>;

const tuples::tagged_tuple_from_typelist<extra_data> initial_extra_data{
    Scalar<DataVector>{{{{123.0}}}}, Scalar<DataVector>{{{{456.0}}}},
    Scalar<double>{{{789.0}}}, "Initial"};

class TestTrigger : public DenseTrigger {
 public:
  TestTrigger() = default;
  explicit TestTrigger(CkMigrateMessage* const msg) : DenseTrigger(msg) {}
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
  WRAPPED_PUPable_decl_template(TestTrigger);  // NOLINT
#pragma GCC diagnostic pop

  // All triggers are evaluated once at the start of the evolution, so
  // we have to handle that call and set up for triggering at the
  // interesting time.
  TestTrigger(const double init_time, const double trigger_time,
              const std::optional<bool>& is_triggered,
              const std::optional<double>& next_trigger)
      : init_time_(init_time),
        trigger_time_(trigger_time),
        is_triggered_(is_triggered),
        next_trigger_(next_trigger) {}

  using is_triggered_return_tags = tmpl::list<>;
  using is_triggered_argument_tags = tmpl::list<Tags::Time>;
  template <typename Metavariables, typename ArrayIndex, typename Component>
  std::optional<bool> is_triggered(
      Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, const Component* /*component*/,
      const double time) const {
    if (time == trigger_time_) {
      return is_triggered_;
    } else {
      // First initialization call.
      CHECK(time == init_time_);
      return false;
    }
  }

  using next_check_time_return_tags = tmpl::list<>;
  using next_check_time_argument_tags = tmpl::list<Tags::Time>;
  template <typename Metavariables, typename ArrayIndex, typename Component>
  std::optional<double> next_check_time(
      Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/, const Component* /*component*/,
      const double time) const {
    if (time == trigger_time_) {
      return next_trigger_;
    } else {
      // First initialization call.
      CHECK(time == init_time_);
      return trigger_time_;
    }
  }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override {
    DenseTrigger::pup(p);
    p | init_time_;
    p | trigger_time_;
    p | is_triggered_;
    p | next_trigger_;
  }

 private:
  double init_time_ = std::numeric_limits<double>::signaling_NaN();
  double trigger_time_ = std::numeric_limits<double>::signaling_NaN();
  std::optional<bool> is_triggered_{};
  std::optional<double> next_trigger_{};
};

PUP::able::PUP_ID TestTrigger::my_PUP_ID = 0;  // NOLINT

struct TestEvent : public Event {
  TestEvent() = default;
  explicit TestEvent(CkMigrateMessage* const /*msg*/) {}
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
  WRAPPED_PUPable_decl_template(TestEvent);  // NOLINT
#pragma GCC diagnostic pop

  explicit TestEvent(const bool needs_evolved_variables)
      : needs_evolved_variables_(needs_evolved_variables) {}

  using compute_tags_for_observation_box = tmpl::list<>;

  using return_tags = tmpl::list<>;
  // Because of a poor choice in the argument order, operator() cannot
  // take a parameter pack of arguments, so we pull out the objects
  // ourselves.
  using argument_tags = tmpl::list<Tags::DataBox>;

  template <typename DbTags, typename Metavariables, typename ArrayIndex,
            typename Component>
  void operator()(const db::DataBox<DbTags>& box,
                  Parallel::GlobalCache<Metavariables>& /*cache*/,
                  const ArrayIndex& /*array_index*/,
                  const Component* const /*meta*/,
                  const ObservationValue& /*observation_value*/) const {
    tmpl::as_pack<all_data>([&](auto... tags_v) {
      calls.emplace_back(db::get<tmpl::type_from<decltype(tags_v)>>(box)...);
    });
  }

  using is_ready_argument_tags = tmpl::list<>;

  template <typename Metavariables, typename ArrayIndex, typename Component>
  bool is_ready(Parallel::GlobalCache<Metavariables>& /*cache*/,
                const ArrayIndex& /*array_index*/,
                const Component* const /*meta*/) const {
    // We use the triggers to control readiness for this test.
    return true;
  }

  bool needs_evolved_variables() const override {
    return needs_evolved_variables_;
  }

  // `modifications` are pairs of Tag{} and functional to compute tag
  // from the evolved variables.
  template <typename... Modifications>
  static void check_calls(
      const std::vector<std::pair<double, EvolvedVariables>>& expected,
      Modifications... modifications) {
    CAPTURE(get_output(expected));
    CAPTURE(get_output(calls));
    REQUIRE(calls.size() == expected.size());
    for (size_t i = 0; i < expected.size(); ++i) {
      CHECK(get<::Tags::Time>(calls[i]) == expected[i].first);
      const auto& expected_evolved = get<EvolvedVar>(expected[i].second);

      const auto modify = [&](auto tag_v, auto expected_value) {
        using tag = decltype(tag_v);
        [[maybe_unused]] const auto apply_modification =
            [&](const auto& modification) {
              if constexpr (std::is_same_v<decltype(modification.first), tag>) {
                expected_value = std::decay_t<decltype(expected_value)>(
                    modification.second(expected_evolved));
              }
              return 0;
            };
        expand_pack(apply_modification(modifications)...);
        return expected_value;
      };

      CHECK(get<EvolvedVar>(calls[i]) ==
            modify(EvolvedVar{}, expected_evolved));
      tmpl::for_each<extra_data>([&](auto tag_v) {
        using tag = tmpl::type_from<decltype(tag_v)>;
        CHECK(get<tag>(calls[i]) ==
              modify(tag{}, get<tag>(initial_extra_data)));
      });
    }
    calls.clear();
  }

 private:
  bool needs_evolved_variables_ = false;

  static std::vector<DataTuple> calls;
};

std::vector<DataTuple> TestEvent::calls{};

PUP::able::PUP_ID TestEvent::my_PUP_ID = 0;  // NOLINT

struct System {
  using variables_tag = Tags::Variables<tmpl::list<EvolvedVar>>;
};

template <typename Metavariables>
struct Component {
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = int;
  using variables_tag = typename metavariables::system::variables_tag;
  using simple_tags_from_options =
      tmpl::push_front<extra_data, Tags::TimeStepId, Tags::TimeStep, Tags::Time,
                       ::Tags::PreviousTriggerTime, variables_tag,
                       Tags::HistoryEvolvedVariables<variables_tag>,
                       ::Tags::EventsAndDenseTriggers,
                       domain::Tags::NeighborMesh<1>, domain::Tags::Element<1>>;
  using compute_tags = time_stepper_ref_tags<TimeStepper>;

  using phase_dependent_action_list = tmpl::list<
      Parallel::PhaseActions<Parallel::Phase::Initialization,
                             tmpl::list<ActionTesting::InitializeDataBox<
                                 tmpl::list<>, compute_tags>>>,
      Parallel::PhaseActions<
          Parallel::Phase::Testing,
          tmpl::list<evolution::Actions::RunEventsAndDenseTriggers<
              typename Metavariables::postprocessors>>>>;
};

template <typename Postprocessors>
struct Metavariables {
  using postprocessors = Postprocessors;
  using system = System;
  using component_list = tmpl::list<Component<Metavariables>>;
  using const_global_cache_tags =
      tmpl::list<Tags::ConcreteTimeStepper<TimeStepper>>;
  struct factory_creation
      : tt::ConformsTo<Options::protocols::FactoryCreation> {
    using factory_classes =
        tmpl::map<tmpl::pair<DenseTrigger, tmpl::list<TestTrigger>>,
                  tmpl::pair<Event, tmpl::list<TestEvent>>>;
  };
};

template <typename Metavariables>
bool run_if_ready(
    const gsl::not_null<ActionTesting::MockRuntimeSystem<Metavariables>*>
        runner) {
  using component = Component<Metavariables>;

  const auto get_data = [&runner]() {
    return tmpl::as_pack<all_data>([&runner](auto... tags_v) {
      return DataTuple(
          ActionTesting::get_databox_tag<component,
                                         tmpl::type_from<decltype(tags_v)>>(
              *runner, 0)...);
    });
  };

  const auto data_before = get_data();
  const bool was_ready =
      ActionTesting::next_action_if_ready<component>(runner, 0);
  const auto data_after = get_data();
  CHECK(data_before == data_after);
  return was_ready;
}

template <typename TestCase>
void test(const bool time_runs_forward) {
  using metavars = typename TestCase::metavariables;
  using MockRuntimeSystem = typename TestCase::MockRuntimeSystem;
  using component = Component<metavars>;
  using system = typename metavars::system;
  using variables_tag = typename system::variables_tag;
  using VarsType = typename variables_tag::type;
  using DtVarsType =
      typename db::add_tag_prefix<::Tags::dt, variables_tag>::type;
  using History = TimeSteppers::History<VarsType>;

  const Slab slab(0.0, 4.0);
  const TimeStepId time_step_id(time_runs_forward, 0,
                                slab.start() + slab.duration() / 2);
  const TimeDelta exact_step_size =
      (time_runs_forward ? 1 : -1) * slab.duration() / 4;
  const double start_time = time_step_id.step_time().value();
  const double step_size = exact_step_size.value();
  const double step_center = start_time + 0.5 * step_size;
  const double done_time = (time_runs_forward ? 1.0 : -1.0) *
                           std::numeric_limits<double>::infinity();
  const VarsType initial_vars{1, 8.0};
  const VarsType stored_vars{1, 123.0};
  const DtVarsType deriv_vars{1, 1.0};
  const VarsType center_vars = initial_vars + 0.5 * step_size * deriv_vars;

  const auto set_up_component =
      [&deriv_vars, &exact_step_size, &initial_vars, &start_time, &stored_vars,
       &time_step_id](
          const gsl::not_null<MockRuntimeSystem*> runner,
          const std::vector<std::tuple<double, std::optional<bool>,
                                       std::optional<double>, bool>>&
              triggers) {
        History history(1);
        history.insert(time_step_id, initial_vars, deriv_vars);

        EventsAndDenseTriggers::ConstructionType events_and_dense_triggers{};
        events_and_dense_triggers.reserve(triggers.size());
        for (auto [trigger_time, is_triggered, next_trigger,
                   needs_evolved_variables] : triggers) {
          events_and_dense_triggers.push_back(
              {std::make_unique<TestTrigger>(start_time, trigger_time,
                                             is_triggered, next_trigger),
               make_vector<std::unique_ptr<Event>>(
                   std::make_unique<TestEvent>(needs_evolved_variables))});
        }

        tmpl::as_pack<extra_data>([&](auto... tags_v) {
          ActionTesting::emplace_array_component_and_initialize<component>(
              runner, ActionTesting::NodeId{0}, ActionTesting::LocalCoreId{0},
              0, {}, time_step_id, exact_step_size, start_time,
              std::optional<double>{}, stored_vars, std::move(history),
              EventsAndDenseTriggers(std::move(events_and_dense_triggers)),
              typename domain::Tags::NeighborMesh<1>::type{},
              Element<1>{ElementId<1>{0}, {}},
              get<tmpl::type_from<decltype(tags_v)>>(initial_extra_data)...);
        });
        ActionTesting::set_phase(runner, Parallel::Phase::Testing);
      };

  // Tests start here

  // Nothing should happen in self-start
  {
    // This isn't a valid time for the trigger to reschedule to (it is
    // in the past), but the triggers should be completely ignored in
    // this check.
    const double invalid_time = start_time - step_size;

    MockRuntimeSystem runner{
        {std::make_unique<TimeSteppers::AdamsBashforth>(1)}};
    set_up_component(&runner,
                     {{invalid_time, std::nullopt, std::nullopt, false}});
    {
      auto& box =
          ActionTesting::get_databox<component>(make_not_null(&runner), 0);
      db::mutate<Tags::TimeStepId>(
          [](const gsl::not_null<TimeStepId*> id) {
            *id = TimeStepId(id->time_runs_forward(), -1, id->step_time());
          },
          make_not_null(&box));
    }
    CHECK(run_if_ready(make_not_null(&runner)));
    TestEvent::check_calls({});
  }

  // No triggers
  {
    MockRuntimeSystem runner{
        {std::make_unique<TimeSteppers::AdamsBashforth>(1)}};
    set_up_component(&runner, {});
    CHECK(run_if_ready(make_not_null(&runner)));
    TestEvent::check_calls({});
  }

  // Triggers too far in the future
  const auto check_not_reached = [&set_up_component, &start_time, &step_size](
                                     const std::optional<bool>& is_triggered) {
    MockRuntimeSystem runner{
        {std::make_unique<TimeSteppers::AdamsBashforth>(1)}};
    set_up_component(&runner, {{start_time + 1.5 * step_size, is_triggered,
                                std::nullopt, false}});
    CHECK(run_if_ready(make_not_null(&runner)));
    TestEvent::check_calls({});
  };
  check_not_reached(std::nullopt);
  check_not_reached(true);
  check_not_reached(false);

  // Triggers at the end of the step (should not run)
  const auto check_step_end = [&set_up_component, &start_time, &step_size](
                                  const std::optional<bool>& is_triggered) {
    MockRuntimeSystem runner{
        {std::make_unique<TimeSteppers::AdamsBashforth>(1)}};
    set_up_component(
        &runner, {{start_time + step_size, is_triggered, std::nullopt, false}});
    CHECK(run_if_ready(make_not_null(&runner)));
    TestEvent::check_calls({});
  };
  check_step_end(std::nullopt);
  check_step_end(true);
  check_step_end(false);

  // Trigger isn't ready
  {
    MockRuntimeSystem runner{
        {std::make_unique<TimeSteppers::AdamsBashforth>(1)}};
    set_up_component(&runner, {{step_center, std::nullopt, start_time, false}});
    CHECK(not run_if_ready(make_not_null(&runner)));
    TestEvent::check_calls({});
  }

  // Variables not needed
  const auto check_not_needed = [&](const bool reschedule) {
    MockRuntimeSystem runner{
        {std::make_unique<TimeSteppers::AdamsBashforth>(1)}};
    const auto next_check =
        reschedule ? std::optional{done_time} : std::nullopt;
    set_up_component(&runner, {{step_center, true, next_check, false}});
    CHECK(run_if_ready(make_not_null(&runner)) == reschedule);
    TestEvent::check_calls({{step_center, stored_vars}});
  };
  check_not_needed(true);
  check_not_needed(false);

  // Variables not needed at initial time
  const auto check_not_needed_initial = [&](const bool reschedule) {
    MockRuntimeSystem runner{
        {std::make_unique<TimeSteppers::AdamsBashforth>(1)}};
    const auto next_check =
        reschedule ? std::optional{done_time} : std::nullopt;
    set_up_component(&runner, {{start_time, true, next_check, false}});
    CHECK(run_if_ready(make_not_null(&runner)) == reschedule);
    TestEvent::check_calls({{start_time, stored_vars}});
  };
  check_not_needed_initial(true);
  check_not_needed_initial(false);

  // Variables needed
  const auto check_needed = [&](const bool reschedule) {
    MockRuntimeSystem runner{
        {std::make_unique<TimeSteppers::AdamsBashforth>(1)}};
    const auto next_check =
        reschedule ? std::optional{done_time} : std::nullopt;
    set_up_component(&runner, {{step_center, true, next_check, true}});
    TestCase::check_dense(&runner, reschedule, {{step_center, center_vars}});
  };
  check_needed(true);
  check_needed(false);

  // Variables needed at initial time
  const auto check_needed_initial = [&](const bool reschedule) {
    MockRuntimeSystem runner{
        {std::make_unique<TimeSteppers::AdamsBashforth>(1)}};
    const auto next_check =
        reschedule ? std::optional{done_time} : std::nullopt;
    set_up_component(&runner, {{start_time, true, next_check, true}});
    TestCase::check_dense(&runner, reschedule, {{start_time, initial_vars}});
  };
  check_needed_initial(true);
  check_needed_initial(false);

  // Missing dense output data
  const auto check_missing_dense_data = [&](const bool data_needed) {
    MockRuntimeSystem runner{
        {std::make_unique<TimeSteppers::Rk3HesthavenSsp>()}};
    set_up_component(&runner, {{step_center, true, done_time, data_needed}});
    {
      auto& box =
          ActionTesting::get_databox<component>(make_not_null(&runner), 0);
      db::mutate<Tags::HistoryEvolvedVariables<variables_tag>>(
          [&deriv_vars, &initial_vars, &time_step_id](
              const gsl::not_null<History*> history,
              const TimeStepper& time_stepper) {
            *history = History(3);
            history->insert(time_step_id, initial_vars, deriv_vars);
            history->insert(
                time_stepper.next_time_id(
                    time_step_id,
                    (time_step_id.time_runs_forward() ? 1 : -1) *
                        time_step_id.step_time().slab().duration() / 4),
                initial_vars, deriv_vars);
          },
          make_not_null(&box), db::get<Tags::TimeStepper<TimeStepper>>(box));
    }
    CHECK(run_if_ready(make_not_null(&runner)));
    if (data_needed) {
      TestEvent::check_calls({});
    } else {
      // If we don't need the data, it shouldn't matter whether it is missing.
      TestEvent::check_calls({{step_center, stored_vars}});
    }
  };
  check_missing_dense_data(true);
  check_missing_dense_data(false);

  // Multiple triggers
  {
    const double second_trigger = start_time + 0.75 * step_size;
    MockRuntimeSystem runner{
        {std::make_unique<TimeSteppers::AdamsBashforth>(1)}};
    set_up_component(&runner, {{step_center, true, done_time, true},
                               {second_trigger, true, done_time, true}});
    TestCase::check_dense(
        &runner, true,
        {{step_center, center_vars},
         {second_trigger, initial_vars + 0.75 * step_size * deriv_vars}});
  }
}

namespace test_postprocessors {
struct SetA {
  using return_tags =
      tmpl::list<PostprocessedVar<Scalar<DataVector>, labels::A>>;
  using argument_tags = tmpl::list<EvolvedVar>;
  static void apply(const gsl::not_null<Scalar<DataVector>*> postprocessed_a,
                    const Scalar<DataVector>& evolved) {
    get(*postprocessed_a) = 2.0 * get(evolved);
  }
};

struct SetAB {
  using return_tags =
      tmpl::list<PostprocessedVar<Scalar<DataVector>, labels::A>,
                 PostprocessedVar<Scalar<DataVector>, labels::B>>;
  using argument_tags = tmpl::list<EvolvedVar>;
  static void apply(const gsl::not_null<Scalar<DataVector>*> postprocessed_a,
                    const gsl::not_null<Scalar<DataVector>*> postprocessed_b,
                    const Scalar<DataVector>& evolved) {
    get(*postprocessed_a) = 3.0 * get(evolved);
    get(*postprocessed_b) = 4.0 * get(evolved);
  }
};

struct SetDouble {
  using return_tags = tmpl::list<PostprocessedVar<Scalar<double>>>;
  using argument_tags = tmpl::list<EvolvedVar>;
  static void apply(const gsl::not_null<Scalar<double>*> postprocessed_double,
                    const Scalar<DataVector>& evolved) {
    get(*postprocessed_double) = 5.0 * get(evolved)[0];
  }

  // Test is_ready
  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ParallelComponent>
  static bool is_ready(
      const gsl::not_null<db::DataBox<DbTagsList>*> /*box*/,
      const gsl::not_null<tuples::TaggedTuple<InboxTags...>*> /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/,
      const ParallelComponent* const /*component*/) {
    return true;
  }
};

struct SetDoubleAndString {
  using return_tags = tmpl::list<PostprocessedVar<Scalar<double>>,
                                 PostprocessedVar<std::string>>;
  using argument_tags = tmpl::list<EvolvedVar>;
  static void apply(const gsl::not_null<Scalar<double>*> postprocessed_double,
                    const gsl::not_null<std::string*> postprocessed_string,
                    const Scalar<DataVector>& evolved) {
    get(*postprocessed_double) = 6.0 * get(evolved)[0];
    *postprocessed_string = "Processed";
  }
};

struct ModifyEvolved {
  using return_tags = tmpl::list<EvolvedVar>;
  using argument_tags = tmpl::list<>;
  static void apply(const gsl::not_null<Scalar<DataVector>*> evolved) {
    get(*evolved) *= -1.0;
  }
};

struct NotReady {
  using return_tags = tmpl::list<>;
  using argument_tags = tmpl::list<>;
  static void apply() {}

  template <typename DbTagsList, typename... InboxTags, typename Metavariables,
            typename ArrayIndex, typename ParallelComponent>
  static bool is_ready(
      const gsl::not_null<db::DataBox<DbTagsList>*> /*box*/,
      const gsl::not_null<tuples::TaggedTuple<InboxTags...>*> /*inboxes*/,
      Parallel::GlobalCache<Metavariables>& /*cache*/,
      const ArrayIndex& /*array_index*/,
      const ParallelComponent* const /*component*/) {
    return false;
  }
};
}  // namespace test_postprocessors

namespace test_cases {
struct NoPostprocessors {
  using postprocessors = tmpl::list<>;
  using metavariables = Metavariables<postprocessors>;
  using MockRuntimeSystem = ActionTesting::MockRuntimeSystem<metavariables>;
  static void check_dense(
      const gsl::not_null<MockRuntimeSystem*> runner, const bool should_run,
      const std::vector<std::pair<double, EvolvedVariables>>& expected_calls) {
    CHECK(run_if_ready(runner) == should_run);
    TestEvent::check_calls(expected_calls);
  }
};

struct NotReady {
  using postprocessors = tmpl::list<test_postprocessors::NotReady>;
  using metavariables = Metavariables<postprocessors>;
  using MockRuntimeSystem = ActionTesting::MockRuntimeSystem<metavariables>;
  static void check_dense(
      const gsl::not_null<MockRuntimeSystem*> runner, const bool /*should_run*/,
      const std::vector<
          std::pair<double, EvolvedVariables>>& /*expected_calls*/) {
    CHECK(not run_if_ready(runner));
    TestEvent::check_calls({});
  }
};

struct PostprocessA {
  using postprocessors =
      tmpl::list<AlwaysReadyPostprocessor<test_postprocessors::SetA>>;
  using metavariables = Metavariables<postprocessors>;
  using MockRuntimeSystem = ActionTesting::MockRuntimeSystem<metavariables>;
  static void check_dense(
      const gsl::not_null<MockRuntimeSystem*> runner, const bool should_run,
      const std::vector<std::pair<double, EvolvedVariables>>& expected_calls) {
    CHECK(run_if_ready(runner) == should_run);
    // We define the lambda out-of-line for nvcc compatibility
    auto multiplier = [](const Scalar<DataVector>& v) { return 2.0 * get(v); };
    TestEvent::check_calls(
        expected_calls,
        std::pair{PostprocessedVar<Scalar<DataVector>, labels::A>{},
                  std::move(multiplier)});
  }
};

struct PostprocessAll {
  // Test setting the same thing multiple times
  using postprocessors = tmpl::list<
      AlwaysReadyPostprocessor<test_postprocessors::SetAB>,
      AlwaysReadyPostprocessor<test_postprocessors::SetA>,
      AlwaysReadyPostprocessor<test_postprocessors::SetDoubleAndString>,
      test_postprocessors::SetDouble>;
  using metavariables = Metavariables<postprocessors>;
  using MockRuntimeSystem = ActionTesting::MockRuntimeSystem<metavariables>;
  static void check_dense(
      const gsl::not_null<MockRuntimeSystem*> runner, const bool should_run,
      const std::vector<std::pair<double, EvolvedVariables>>& expected_calls) {
    CHECK(run_if_ready(runner) == should_run);
    // We define the lambdas out-of-line for nvcc compatibility
    auto multiplier0 = [](const Scalar<DataVector>& v) { return 2.0 * get(v); };
    auto multiplier1 = [](const Scalar<DataVector>& v) { return 4.0 * get(v); };
    auto multiplier2 = [](const Scalar<DataVector>& v) {
      return 5.0 * get(v)[0];
    };
    auto multiplier3 = [](const Scalar<DataVector>& /*v*/) {
      return "Processed";
    };
    TestEvent::check_calls(
        expected_calls,
        std::pair{PostprocessedVar<Scalar<DataVector>, labels::A>{},
                  std::move(multiplier0)},
        std::pair{PostprocessedVar<Scalar<DataVector>, labels::B>{},
                  std::move(multiplier1)},
        std::pair{PostprocessedVar<Scalar<double>>{}, std::move(multiplier2)},
        std::pair{PostprocessedVar<std::string>{}, std::move(multiplier3)});
  }
};

struct PostprocessEvolved {
  using postprocessors =
      tmpl::list<AlwaysReadyPostprocessor<test_postprocessors::ModifyEvolved>,
                 AlwaysReadyPostprocessor<test_postprocessors::SetA>>;
  using metavariables = Metavariables<postprocessors>;
  using MockRuntimeSystem = ActionTesting::MockRuntimeSystem<metavariables>;
  static void check_dense(
      const gsl::not_null<MockRuntimeSystem*> runner, const bool should_run,
      const std::vector<std::pair<double, EvolvedVariables>>& expected_calls) {
    CHECK(run_if_ready(runner) == should_run);
    // We define the lambdas out-of-line for nvcc compatibility
    auto helper0 = [](const Scalar<DataVector>& v) { return -get(v); };
    auto helper1 = [](const Scalar<DataVector>& v) { return -2.0 * get(v); };
    TestEvent::check_calls(
        expected_calls, std::pair{EvolvedVar{}, std::move(helper0)},
        std::pair{PostprocessedVar<Scalar<DataVector>, labels::A>{},
                  std::move(helper1)});
  }
};
}  // namespace test_cases

EventsAndDenseTriggers make_events_and_dense_triggers() {
  const Slab slab(0.0, 4.0);
  const TimeStepId time_step_id(true, 0, slab.start() + slab.duration() / 2);
  const TimeDelta exact_step_size = slab.duration() / 4;
  const double start_time = time_step_id.step_time().value();
  const double step_size = exact_step_size.value();
  const double step_center = start_time + 0.5 * step_size;
  const double second_trigger = start_time + 0.75 * step_size;
  const double done_time = std::numeric_limits<double>::infinity();

  EventsAndDenseTriggers::ConstructionType events_and_dense_triggers{};
  const std::vector<
      std::tuple<double, std::optional<bool>, std::optional<double>, bool>>
      triggers{{step_center, true, done_time, true},
               {second_trigger, true, done_time, true}};
  events_and_dense_triggers.reserve(triggers.size());
  for (auto [trigger_time, is_triggered, next_trigger,
             needs_evolved_variables] : triggers) {
    events_and_dense_triggers.push_back(
        {std::make_unique<TestTrigger>(start_time, trigger_time, is_triggered,
                                       next_trigger),
         make_vector<std::unique_ptr<Event>>(
             std::make_unique<TestEvent>(needs_evolved_variables))});
  }
  return EventsAndDenseTriggers(std::move(events_and_dense_triggers));
}

void test_p_refine() {
  auto box = db::create<db::AddSimpleTags<::Tags::EventsAndDenseTriggers,
                                          ::Tags::PreviousTriggerTime>>(
      make_events_and_dense_triggers(), std::optional<double>{});

  const Element<1> element{ElementId<1>{0}, DirectionMap<1, Neighbors<1>>{}};
  const Mesh<1> mesh{2, Spectral::Basis::Legendre,
                     Spectral::Quadrature::GaussLobatto};
  db::mutate_apply<evolution::Actions::ProjectRunEventsAndDenseTriggers>(
      make_not_null(&box), std::make_pair(mesh, element));

  // There is no comparison operator for EventsAndDenseTriggers
  // const auto expected_events_and_dense_triggers =
  //     make_events_and_dense_triggers();
  // CHECK(db::get<::Tags::EventsAndDenseTriggers>(box) ==
  //       expected_events_and_dense_triggers);
  CHECK(db::get<::Tags::PreviousTriggerTime>(box) == std::nullopt);
}
}  // namespace

SPECTRE_TEST_CASE("Unit.Evolution.RunEventsAndDenseTriggers",
                  "[Unit][Evolution][Actions]") {
  register_classes_with_charm<TimeSteppers::AdamsBashforth,
                              TimeSteppers::Rk3HesthavenSsp>();
  register_factory_classes_with_charm<Metavariables<tmpl::list<>>>();

  for (const auto time_runs_forward : {true, false}) {
    test<test_cases::NoPostprocessors>(time_runs_forward);
    test<test_cases::NotReady>(time_runs_forward);
    test<test_cases::PostprocessA>(time_runs_forward);
    test<test_cases::PostprocessAll>(time_runs_forward);
    test<test_cases::PostprocessEvolved>(time_runs_forward);
  }
  static_assert(tt::assert_conforms_to_v<
                evolution::Actions::ProjectRunEventsAndDenseTriggers,
                amr::protocols::Projector>);
  test_p_refine();
}
