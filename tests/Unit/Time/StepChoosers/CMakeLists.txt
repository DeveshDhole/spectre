# Distributed under the MIT License.
# See LICENSE.txt for details.

set(LIBRARY_SOURCES
  ${LIBRARY_SOURCES}
  StepChoosers/Test_ByBlock.cpp
  StepChoosers/Test_Cfl.cpp
  StepChoosers/Test_Constant.cpp
  StepChoosers/Test_ElementSizeCfl.cpp
  StepChoosers/Test_ErrorControl.cpp
  StepChoosers/Test_FixedLtsRatio.cpp
  StepChoosers/Test_LimitIncrease.cpp
  StepChoosers/Test_Maximum.cpp
  StepChoosers/Test_PreventRapidIncrease.cpp
  StepChoosers/Test_Random.cpp
  StepChoosers/Test_StepToTimes.cpp
  PARENT_SCOPE)
