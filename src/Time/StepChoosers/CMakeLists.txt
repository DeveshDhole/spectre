# Distributed under the MIT License.
# See LICENSE.txt for details.

spectre_target_sources(
  ${LIBRARY}
  PRIVATE
  ByBlock.cpp
  Constant.cpp
  FixedLtsRatio.cpp
  LimitIncrease.cpp
  Maximum.cpp
  PreventRapidIncrease.cpp
  Random.cpp
  StepToTimes.cpp
  )

spectre_target_headers(
  ${LIBRARY}
  INCLUDE_DIRECTORY ${CMAKE_SOURCE_DIR}/src
  HEADERS
  ByBlock.hpp
  Cfl.hpp
  Constant.hpp
  ElementSizeCfl.hpp
  ErrorControl.hpp
  Factory.hpp
  FixedLtsRatio.hpp
  LimitIncrease.hpp
  Maximum.hpp
  PreventRapidIncrease.hpp
  Random.hpp
  StepChooser.hpp
  StepToTimes.hpp
  )
