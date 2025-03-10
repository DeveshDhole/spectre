# Distributed under the MIT License.
# See LICENSE.txt for details.

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  AddTemporalIdsToInterpolationTarget.hpp
  CleanUpInterpolator.hpp
  DumpInterpolatorVolumeData.hpp
  ElementInitInterpPoints.hpp
  ElementReceiveInterpPoints.hpp
  InitializeInterpolationTarget.hpp
  InitializeInterpolator.hpp
  InterpolationTargetReceiveVars.hpp
  InterpolationTargetSendPoints.hpp
  InterpolationTargetVarsFromElement.hpp
  InterpolatorReceivePoints.hpp
  InterpolatorReceiveVolumeData.hpp
  InterpolatorRegisterElement.hpp
  PrintInterpolationTargetForDeadlock.hpp
  PrintInterpolatorForDeadlock.hpp
  SendPointsToInterpolator.hpp
  TryToInterpolate.hpp
  VerifyTemporalIdsAndSendPoints.hpp
  )
