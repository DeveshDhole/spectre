// Distributed under the MIT License.
// See LICENSE.txt for details.

#include "Framework/TestingFramework.hpp"

#include <cstddef>
#include <string>
#include <unordered_map>

#include "Framework/ActionTesting.hpp"
#include "Parallel/Phase.hpp"
#include "Parallel/PhaseDependentActionList.hpp"
#include "ParallelAlgorithms/Interpolation/Actions/InitializeInterpolator.hpp"
#include "ParallelAlgorithms/Interpolation/InterpolatedVars.hpp"
#include "ParallelAlgorithms/Interpolation/Tags.hpp"
#include "PointwiseFunctions/GeneralRelativity/Tags.hpp"
#include "Time/Tags/Time.hpp"
#include "Time/Tags/TimeStepId.hpp"
#include "Utilities/Gsl.hpp"
#include "Utilities/TMPL.hpp"
#include "Utilities/TaggedTuple.hpp"

class DataVector;

namespace {

template <typename Metavariables>
struct mock_interpolator {
  using metavariables = Metavariables;
  using chare_type = ActionTesting::MockArrayChare;
  using array_index = size_t;
  using phase_dependent_action_list = tmpl::list<Parallel::PhaseActions<
      Parallel::Phase::Initialization,
      tmpl::list<intrp::Actions::InitializeInterpolator<
          Metavariables::volume_dim,
          tmpl::list<
              intrp::Tags::VolumeVarsInfo<Metavariables, ::Tags::Time>,
              intrp::Tags::VolumeVarsInfo<Metavariables, ::Tags::TimeStepId>>,
          intrp::Tags::InterpolatedVarsHolders<Metavariables>>>>>;
};

struct Metavariables {
  struct InterpolatorTargetA {
    using temporal_id = ::Tags::TimeStepId;
    using vars_to_interpolate_to_target =
        tmpl::list<gr::Tags::Lapse<DataVector>>;
  };
  static constexpr size_t volume_dim = 3;
  using interpolator_source_vars = tmpl::list<gr::Tags::Lapse<DataVector>>;
  using interpolation_target_tags = tmpl::list<InterpolatorTargetA>;

  using component_list = tmpl::list<mock_interpolator<Metavariables>>;
};

SPECTRE_TEST_CASE("Unit.NumericalAlgorithms.Interpolator.Initialize",
                  "[Unit]") {
  using metavars = Metavariables;
  using component = mock_interpolator<metavars>;
  ActionTesting::MockRuntimeSystem<metavars> runner{{}};
  ActionTesting::set_phase(make_not_null(&runner),
                           Parallel::Phase::Initialization);
  ActionTesting::emplace_component<component>(&runner, 0);
  for (size_t i = 0; i < 2; ++i) {
    ActionTesting::next_action<component>(make_not_null(&runner), 0);
  }
  ActionTesting::set_phase(make_not_null(&runner), Parallel::Phase::Testing);

  CHECK(ActionTesting::get_databox_tag<
            component, ::intrp::Tags::NumberOfElements<metavars::volume_dim>>(
            runner, 0)
            .empty());
  CHECK(ActionTesting::get_databox_tag<
            component, ::intrp::Tags::VolumeVarsInfo<metavars, ::Tags::Time>>(
            runner, 0)
            .empty());
  CHECK(ActionTesting::get_databox_tag<
            component,
            ::intrp::Tags::VolumeVarsInfo<metavars, ::Tags::TimeStepId>>(runner,
                                                                         0)
            .empty());

  const auto& holders = ActionTesting::get_databox_tag<
      component, ::intrp::Tags::InterpolatedVarsHolders<metavars>>(runner, 0);
  // Check that 'holders' has a tag corresponding to
  // 'metavars::InterpolatorTargetA'
  const auto& holder =
      get<intrp::Vars::HolderTag<metavars::InterpolatorTargetA, metavars>>(
          holders);
  CHECK(holder.infos.empty());
  // Check that 'holders' has only one tag.
  CHECK(holders.size() == 1);
}
}  // namespace
