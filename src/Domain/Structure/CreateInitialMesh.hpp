// Distributed under the MIT License.
// See LICENSE.txt for details.

#pragma once

#include <array>
#include <cstddef>
#include <cstdint>
#include <vector>

/// \cond
template <size_t Dim>
struct ElementId;
template <size_t Dim>
class Mesh;
template <size_t Dim>
struct OrientationMap;
namespace Spectral {
enum class Quadrature : uint8_t;
}  // namespace Spectral
/// \endcond

namespace domain::Initialization {
/// \ingroup InitializationGroup
/// \brief Construct the initial Mesh of an Element.
///
/// \details When constructing the Mesh of an Element, pass its id, and use the
/// default argument for orientation.  When constructing the mesh of a
/// neighboring Element (when constructing mortars), pass the id and orientation
/// of the neighbor.
///
/// \param initial_extents the initial extents of each Block in the Domain
/// \param element_id id of an Element or its neighbor
/// \param quadrature the quadrature rule/grid point distribution
/// \param orientation OrientationMap of (neighboring) `element_id`
template <size_t Dim>
Mesh<Dim> create_initial_mesh(
    const std::vector<std::array<size_t, Dim>>& initial_extents,
    const ElementId<Dim>& element_id, Spectral::Quadrature quadrature,
    const OrientationMap<Dim>& orientation =
        OrientationMap<Dim>::create_aligned());
}  // namespace domain::Initialization
