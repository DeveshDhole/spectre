# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY EllipticSystems)

add_spectre_library(${LIBRARY} INTERFACE)

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  GetFluxesComputer.hpp
  GetSourcesComputer.hpp
  )

target_link_libraries(
  ${LIBRARY}
  INTERFACE
  Utilities
  )

add_subdirectory(Elasticity)
add_subdirectory(BnsInitialData)
add_subdirectory(Poisson)
add_subdirectory(Punctures)
add_subdirectory(ScalarGaussBonnet)
add_subdirectory(Xcts)

