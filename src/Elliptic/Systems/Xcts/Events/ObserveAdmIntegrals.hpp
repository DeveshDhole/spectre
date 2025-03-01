// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <optional>
#include <pup.h>
#include <string>
#include <vector>

#include "DataStructures/DataVector.hpp"
#include "DataStructures/Tensor/Tensor.hpp"
#include "Domain/Structure/ElementId.hpp"
#include "Domain/Tags.hpp"
#include "Domain/Tags/FaceNormal.hpp"
#include "Domain/Tags/Faces.hpp"
#include "Elliptic/Systems/Xcts/Tags.hpp"
#include "IO/Observer/GetSectionObservationKey.hpp"
#include "IO/Observer/ObservationId.hpp"
#include "IO/Observer/ObserverComponent.hpp"
#include "IO/Observer/ReductionActions.hpp"
#include "IO/Observer/TypeOfObservation.hpp"
#include "NumericalAlgorithms/Spectral/Mesh.hpp"
#include "Options/String.hpp"
#include "Parallel/ArrayIndex.hpp"
#include "Parallel/GlobalCache.hpp"
#include "Parallel/Reduction.hpp"
#include "Parallel/TypeTraits.hpp"
#include "ParallelAlgorithms/Events/Tags.hpp"
#include "ParallelAlgorithms/EventsAndTriggers/Event.hpp"
#include "Utilities/Serialization/CharmPupable.hpp"

namespace Events {

/// @{
/*!
 * \brief Computes the ADM integrals locally (within one element).
 *
 * To get the total ADM integrals, the results need to be summed over in a
 * reduction.
 */
void local_adm_integrals(
    gsl::not_null<Scalar<double>*> adm_mass,
    gsl::not_null<tnsr::I<double, 3>*> adm_linear_momentum,
    gsl::not_null<Scalar<double>*> adm_angular_momentum_z,
    gsl::not_null<tnsr::I<double, 3>*> center_of_mass,
    const Scalar<DataVector>& conformal_factor,
    const tnsr::i<DataVector, 3>& deriv_conformal_factor,
    const tnsr::ii<DataVector, 3>& conformal_metric,
    const tnsr::II<DataVector, 3>& inv_conformal_metric,
    const tnsr::Ijj<DataVector, 3>& conformal_christoffel_second_kind,
    const tnsr::i<DataVector, 3>& conformal_christoffel_contracted,
    const tnsr::ii<DataVector, 3>& spatial_metric,
    const tnsr::II<DataVector, 3>& inv_spatial_metric,
    const tnsr::ii<DataVector, 3>& extrinsic_curvature,
    const Scalar<DataVector>& trace_extrinsic_curvature,
    const tnsr::I<DataVector, 3, Frame::Inertial>& inertial_coords,
    const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                          Frame::Inertial>& inv_jacobian,
    const Mesh<3>& mesh, const Element<3>& element,
    const DirectionMap<3, tnsr::i<DataVector, 3>>& conformal_face_normals);
/// @}

/// @{
/*!
 * \brief Observe ADM integrals after the XCTS solve.
 *
 * The surface integrals are taken over the outer boundary, which is defined as
 * the domain boundary in the upper logical zeta direction.
 *
 * Writes reduction quantities:
 * - Number of points in the domain
 * - ADM mass
 * - ADM linear momentum
 * - ADM angular momentum (z-component)
 * - Center of mass
 */
template <typename ArraySectionIdTag = void>
class ObserveAdmIntegrals : public Event {
 private:
  using ReductionData = Parallel::ReductionData<
      // Number of points
      Parallel::ReductionDatum<size_t, funcl::Plus<>>,
      // ADM Mass
      Parallel::ReductionDatum<double, funcl::Plus<>>,
      // ADM Linear Momentum (x-component)
      Parallel::ReductionDatum<double, funcl::Plus<>>,
      // ADM Linear Momentum (y-component)
      Parallel::ReductionDatum<double, funcl::Plus<>>,
      // ADM Linear Momentum (z-component)
      Parallel::ReductionDatum<double, funcl::Plus<>>,
      // ADM Angular Momentum (z-component)
      Parallel::ReductionDatum<double, funcl::Plus<>>,
      // Center of Mass (x-component)
      Parallel::ReductionDatum<double, funcl::Plus<>, funcl::Divides<>,
                               std::index_sequence<1>>,
      // Center of Mass (y-component)
      Parallel::ReductionDatum<double, funcl::Plus<>, funcl::Divides<>,
                               std::index_sequence<1>>,
      // Center of Mass (z-component)
      Parallel::ReductionDatum<double, funcl::Plus<>, funcl::Divides<>,
                               std::index_sequence<1>>>;

 public:
  /// \cond
  explicit ObserveAdmIntegrals(CkMigrateMessage* msg) : Event(msg) {}
  using PUP::able::register_constructor;
  WRAPPED_PUPable_decl_template(ObserveAdmIntegrals);  // NOLINT
  /// \endcond

  using options = tmpl::list<>;
  static constexpr Options::String help =
      "Observe ADM integrals after the XCTS solve.\n"
      "\n"
      "Writes reduction quantities:\n"
      "- Number of points in the domain\n"
      "- ADM mass\n"
      "- ADM linear momentum\n"
      "- ADM angular momentum (z-component)\n"
      "- Center of mass";

  ObserveAdmIntegrals() = default;

  using observed_reduction_data_tags =
      observers::make_reduction_data_tags<tmpl::list<ReductionData>>;

  using compute_tags_for_observation_box = tmpl::list<>;

  using return_tags = tmpl::list<>;

  using argument_tags = tmpl::list<
      Xcts::Tags::ConformalFactor<DataVector>,
      ::Tags::deriv<Xcts::Tags::ConformalFactor<DataVector>, tmpl::size_t<3>,
                    Frame::Inertial>,
      Xcts::Tags::ConformalMetric<DataVector, 3, Frame::Inertial>,
      Xcts::Tags::InverseConformalMetric<DataVector, 3, Frame::Inertial>,
      Xcts::Tags::ConformalChristoffelSecondKind<DataVector, 3,
                                                 Frame::Inertial>,
      Xcts::Tags::ConformalChristoffelContracted<DataVector, 3,
                                                 Frame::Inertial>,
      gr::Tags::SpatialMetric<DataVector, 3, Frame::Inertial>,
      gr::Tags::InverseSpatialMetric<DataVector, 3, Frame::Inertial>,
      gr::Tags::ExtrinsicCurvature<DataVector, 3, Frame::Inertial>,
      gr::Tags::TraceExtrinsicCurvature<DataVector>,
      domain::Tags::Coordinates<3, Frame::Inertial>,
      domain::Tags::InverseJacobian<3, Frame::ElementLogical, Frame::Inertial>,
      domain::Tags::Mesh<3>, domain::Tags::Element<3>,
      domain::Tags::Faces<3, domain::Tags::FaceNormal<3>>,
      ::Tags::ObservationBox>;

  template <typename DataBoxType, typename ComputeTagsList,
            typename Metavariables, typename ArrayIndex,
            typename ParallelComponent>
  void operator()(
      const Scalar<DataVector>& conformal_factor,
      const tnsr::i<DataVector, 3>& deriv_conformal_factor,
      const tnsr::ii<DataVector, 3>& conformal_metric,
      const tnsr::II<DataVector, 3>& inv_conformal_metric,
      const tnsr::Ijj<DataVector, 3>& conformal_christoffel_second_kind,
      const tnsr::i<DataVector, 3>& conformal_christoffel_contracted,
      const tnsr::ii<DataVector, 3>& spatial_metric,
      const tnsr::II<DataVector, 3>& inv_spatial_metric,
      const tnsr::ii<DataVector, 3>& extrinsic_curvature,
      const Scalar<DataVector>& trace_extrinsic_curvature,
      const tnsr::I<DataVector, 3, Frame::Inertial>& inertial_coords,
      const InverseJacobian<DataVector, 3, Frame::ElementLogical,
                            Frame::Inertial>& inv_jacobian,
      const Mesh<3>& mesh, const Element<3>& element,
      const DirectionMap<3, tnsr::i<DataVector, 3>>& conformal_face_normals,
      const ObservationBox<DataBoxType, ComputeTagsList>& box,
      Parallel::GlobalCache<Metavariables>& cache,
      const ArrayIndex& array_index, const ParallelComponent* const /*meta*/,
      const ObservationValue& observation_value) const {
    // Skip observation on elements that are not part of a section
    const std::optional<std::string> section_observation_key =
        observers::get_section_observation_key<ArraySectionIdTag>(box);
    if (not section_observation_key.has_value()) {
      return;
    }
    const std::string subfile_path = subfile_path_ + *section_observation_key;

    Scalar<double> adm_mass;
    tnsr::I<double, 3> adm_linear_momentum;
    Scalar<double> adm_angular_momentum_z;
    tnsr::I<double, 3> center_of_mass;
    local_adm_integrals(
        make_not_null(&adm_mass), make_not_null(&adm_linear_momentum),
        make_not_null(&adm_angular_momentum_z), make_not_null(&center_of_mass),
        conformal_factor, deriv_conformal_factor, conformal_metric,
        inv_conformal_metric, conformal_christoffel_second_kind,
        conformal_christoffel_contracted, spatial_metric, inv_spatial_metric,
        extrinsic_curvature, trace_extrinsic_curvature, inertial_coords,
        inv_jacobian, mesh, element, conformal_face_normals);

    // Save components of linear momentum as reduction data
    ReductionData reduction_data{
        mesh.number_of_grid_points(), get(adm_mass),
        get<0>(adm_linear_momentum),  get<1>(adm_linear_momentum),
        get<2>(adm_linear_momentum),  get(adm_angular_momentum_z),
        get<0>(center_of_mass),       get<1>(center_of_mass),
        get<2>(center_of_mass)};
    std::vector<std::string> legend{
        "NumberOfPoints",      "AdmMass",
        "AdmLinearMomentum_x", "AdmLinearMomentum_y",
        "AdmLinearMomentum_z", "AdmAngularMomentum_z",
        "CenterOfMass_x",      "CenterOfMass_y",
        "CenterOfMass_z"};

    // Get information required for reduction
    auto& local_observer = *Parallel::local_branch(
        Parallel::get_parallel_component<
            tmpl::conditional_t<Parallel::is_nodegroup_v<ParallelComponent>,
                                observers::ObserverWriter<Metavariables>,
                                observers::Observer<Metavariables>>>(cache));
    observers::ObservationId observation_id{observation_value.value,
                                            subfile_path + ".dat"};
    Parallel::ArrayComponentId array_component_id{
        std::add_pointer_t<ParallelComponent>{nullptr},
        Parallel::ArrayIndex<ElementId<3>>(array_index)};

    // Send reduction action
    if constexpr (Parallel::is_nodegroup_v<ParallelComponent>) {
      Parallel::threaded_action<
          observers::ThreadedActions::CollectReductionDataOnNode>(
          local_observer, std::move(observation_id),
          std::move(array_component_id), subfile_path, std::move(legend),
          std::move(reduction_data));
    } else {
      Parallel::simple_action<observers::Actions::ContributeReductionData>(
          local_observer, std::move(observation_id),
          std::move(array_component_id), subfile_path, std::move(legend),
          std::move(reduction_data));
    }
  }

  using observation_registration_tags = tmpl::list<::Tags::DataBox>;

  template <typename DbTagsList>
  std::optional<
      std::pair<observers::TypeOfObservation, observers::ObservationKey>>
  get_observation_type_and_key_for_registration(
      const db::DataBox<DbTagsList>& box) const {
    const std::optional<std::string> section_observation_key =
        observers::get_section_observation_key<ArraySectionIdTag>(box);
    if (not section_observation_key.has_value()) {
      return std::nullopt;
    }
    return {{observers::TypeOfObservation::Reduction,
             observers::ObservationKey(
                 subfile_path_ + section_observation_key.value() + ".dat")}};
  }

  using is_ready_argument_tags = tmpl::list<>;

  template <typename Metavariables, typename ArrayIndex, typename Component>
  bool is_ready(Parallel::GlobalCache<Metavariables>& /*cache*/,
                const ArrayIndex& /*array_index*/,
                const Component* const /*meta*/) const {
    return true;
  }

  bool needs_evolved_variables() const override { return false; }

  // NOLINTNEXTLINE(google-runtime-references)
  void pup(PUP::er& p) override {
    Event::pup(p);
    p | subfile_path_;
  }

 private:
  std::string subfile_path_ = "/AdmIntegrals";
};
/// @}

/// \cond
template <typename ArraySectionIdTag>
PUP::able::PUP_ID ObserveAdmIntegrals<ArraySectionIdTag>::my_PUP_ID =
    0;  // NOLINT
/// \endcond

}  // namespace Events
