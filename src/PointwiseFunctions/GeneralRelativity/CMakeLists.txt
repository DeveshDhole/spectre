# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY GeneralRelativity)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  Christoffel.cpp
  DerivativeSpatialMetric.cpp
  DerivativesOfSpacetimeMetric.cpp
  ExtrinsicCurvature.cpp
  GeodesicAcceleration.cpp
  InterfaceNullNormal.cpp
  InverseSpacetimeMetric.cpp
  KerrSchildCoords.cpp
  Lapse.cpp
  ProjectionOperators.cpp
  Psi4.cpp
  Psi4Real.cpp
  Ricci.cpp
  Shift.cpp
  SpacetimeDerivativeOfGothG.cpp
  SpacetimeMetric.cpp
  SpacetimeNormalOneForm.cpp
  SpacetimeNormalVector.cpp
  SpatialMetric.cpp
  TimeDerivativeOfSpacetimeMetric.cpp
  TimeDerivativeOfSpatialMetric.cpp
  TortoiseCoordinates.cpp
  WeylElectric.cpp
  WeylMagnetic.cpp
  WeylPropagating.cpp
  WeylTypeD1.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  Christoffel.hpp
  DerivativeSpatialMetric.hpp
  DerivativesOfSpacetimeMetric.hpp
  DetAndInverseSpatialMetric.hpp
  ExtrinsicCurvature.hpp
  GeodesicAcceleration.hpp
  InterfaceNullNormal.hpp
  InverseSpacetimeMetric.hpp
  KerrSchildCoords.hpp
  Lapse.hpp
  ProjectionOperators.hpp
  Ricci.hpp
  Psi4.hpp
  Psi4Real.hpp
  Shift.hpp
  SpacetimeDerivativeOfGothG.hpp
  SpacetimeMetric.hpp
  SpacetimeNormalOneForm.hpp
  SpacetimeNormalVector.hpp
  SpatialMetric.hpp
  Tags.hpp
  TagsDeclarations.hpp
  TimeDerivativeOfSpacetimeMetric.hpp
  TimeDerivativeOfSpatialMetric.hpp
  TortoiseCoordinates.hpp
  WeylElectric.hpp
  WeylMagnetic.hpp
  WeylPropagating.hpp
  WeylTypeD1.hpp
  )

target_link_libraries(
  ${LIBRARY}
  PUBLIC
  DataStructures
  Domain
  Utilities
  PRIVATE
  RootFinding
  INTERFACE
  ErrorHandling
  )

add_subdirectory(GeneralizedHarmonic)
add_subdirectory(Python)
add_subdirectory(Surfaces)
add_subdirectory(Tags)
