# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY LinearOperators)

add_spectre_library(${LIBRARY})

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  CoefficientTransforms.cpp
  DefiniteIntegral.cpp
  Divergence.cpp
  ExponentialFilter.cpp
  IndefiniteIntegral.cpp
  Linearize.cpp
  PartialDerivatives.cpp
  MeanValue.cpp
  PowerMonitors.cpp
  WeakDivergence.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  CoefficientTransforms.hpp
  DefiniteIntegral.hpp
  Divergence.hpp
  Divergence.tpp
  ExponentialFilter.hpp
  IndefiniteIntegral.hpp
  Linearize.hpp
  MeanValue.hpp
  PartialDerivatives.hpp
  PartialDerivatives.tpp
  PowerMonitors.hpp
  WeakDivergence.hpp
  )

target_link_libraries(
  ${LIBRARY}
  PUBLIC
  BLAS::BLAS
  DataStructures
  ErrorHandling
  Options
  Serialization
  Spectral
  SphericalHarmonics
  Utilities
  INTERFACE
  Domain
  DomainStructure
  )

add_subdirectory(Python)
