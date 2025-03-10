# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY "Test_GhValenciaDivClean")

set(LIBRARY_SOURCES
  Actions/Test_SetInitialData.cpp
  BoundaryConditions/Test_ConstraintPreservingFreeOutflow.cpp
  BoundaryConditions/Test_DirichletAnalytic.cpp
  BoundaryConditions/Test_DirichletFreeOutflow.cpp
  BoundaryCorrections/Test_ProductOfCorrections.cpp
  FiniteDifference/Test_BoundaryConditionGhostData.cpp
  FiniteDifference/Test_Derivatives.cpp
  FiniteDifference/Test_FilterOptions.cpp
  FiniteDifference/Test_Filters.cpp
  FiniteDifference/Test_MonotonisedCentral.cpp
  FiniteDifference/Test_PositivityPreservingAdaptiveOrder.cpp
  FiniteDifference/Test_Wcns5z.cpp
  Subcell/Test_FixConservativesAndComputePrims.cpp
  Subcell/Test_NeighborPackagedData.cpp
  Subcell/Test_PrimitiveGhostData.cpp
  Subcell/Test_PrimsAfterRollback.cpp
  Subcell/Test_ResizeAndComputePrimitives.cpp
  Test_Characteristics.cpp
  Test_Constraints.cpp
  Test_SetPiAndPhiFromConstraints.cpp
  Test_StressEnergy.cpp
  Test_Tags.cpp
  Test_TimeDerivativeTerms.cpp
  )

# The following test is disabled for nvcc because nvcc can't handle it
# (it runs out of memory).
if(NOT ${KOKKOS_CXX_COMPILER_ID} STREQUAL "NVIDIA")
  list(APPEND LIBRARY_SOURCES Subcell/Test_TimeDerivative.cpp)
endif()

add_test_library(${LIBRARY} "${LIBRARY_SOURCES}")

target_link_libraries(
  ${LIBRARY}
  PRIVATE
  DomainBoundaryConditionsHelpers
  DomainHelpers
  Framework
  GeneralRelativityHelpers
  GeneralRelativitySolutions
  GeneralizedHarmonic
  GhGrMhdSolutions
  GhValenciaDivClean
  HydroHelpers
  ValenciaDivClean
  )
