# Distributed under the MIT License.
# See LICENSE.txt for details.

spectre_define_test_timeout_factor_option(INPUT_FILE "input file")

option(SPECTRE_INPUT_FILE_TEST_MIN_PRIORITY "Minimum priority of input file \
tests to run. Possible values are: low (not usually run on CI), normal \
(run at least once on CI), high (run always on CI)." "normal")

# Convert priority string to a number.
if (NOT SPECTRE_INPUT_FILE_TEST_MIN_PRIORITY)
  set(SPECTRE_INPUT_FILE_TEST_MIN_PRIORITY "normal")
endif()
string(TOLOWER "${SPECTRE_INPUT_FILE_TEST_MIN_PRIORITY}"
  SPECTRE_INPUT_FILE_TEST_MIN_PRIORITY)
if (${SPECTRE_INPUT_FILE_TEST_MIN_PRIORITY} STREQUAL "high")
  set(SPECTRE_INPUT_FILE_TEST_MIN_PRIORITY 1)
elseif (${SPECTRE_INPUT_FILE_TEST_MIN_PRIORITY} STREQUAL "normal")
  set(SPECTRE_INPUT_FILE_TEST_MIN_PRIORITY 0)
elseif (${SPECTRE_INPUT_FILE_TEST_MIN_PRIORITY} STREQUAL "low")
  set(SPECTRE_INPUT_FILE_TEST_MIN_PRIORITY -1)
else()
  message(FATAL_ERROR "Unknown priority in option "
    "SPECTRE_INPUT_FILE_TEST_MIN_PRIORITY: "
    "${SPECTRE_INPUT_FILE_TEST_MIN_PRIORITY}. "
    "Possible values are: high, normal, low")
endif()

# Environment variables for test
set(_INPUT_FILE_TEST_ENV_VARS "")
# - Disable ASAN's leak sanitizer because Charm++ has false positives
list(APPEND _INPUT_FILE_TEST_ENV_VARS "ASAN_OPTIONS=detect_leaks=0")
# - Set PYTHONPATH to find Python modules
list(APPEND _INPUT_FILE_TEST_ENV_VARS "PYTHONPATH=${PYTHONPATH}")

function(add_single_input_file_test INPUT_FILE EXECUTABLE COMMAND_LINE_ARGS
                                    CHECK_TYPE TIMEOUT EXPECTED_EXIT_CODE)
  # Extract just the name of the input file
  get_filename_component(INPUT_FILE_NAME "${INPUT_FILE}" NAME)

  # Extract the main subdirectory name
  string(FIND "${INPUT_FILE}" "tests/InputFiles/" POSITION_OF_INPUT_FILE_DIR)
  math(EXPR
    POSITION_OF_INPUT_FILE_DIR
    ${POSITION_OF_INPUT_FILE_DIR}+17
    # 17 is the length of "tests/InputFiles/"
    )
  string(SUBSTRING "${INPUT_FILE}" ${POSITION_OF_INPUT_FILE_DIR}
    -1 TEMP)
  string(FIND "${TEMP}" "/" POSITION_OF_SLASH)
  string(SUBSTRING "${TEMP}" 0 ${POSITION_OF_SLASH}
      EXECUTABLE_DIR_NAME)

  # Set tags for the test
  set(TAGS "InputFiles;${EXECUTABLE_DIR_NAME};${CHECK_TYPE}")
  string(TOLOWER "${TAGS}" TAGS)
  # Add the executable name as label, without converting to lower case. This
  # allows running all input file tests for a particular executable.
  list(APPEND TAGS ${EXECUTABLE})

  set(
    CTEST_NAME
    "InputFiles.${EXECUTABLE_DIR_NAME}.${INPUT_FILE_NAME}.${CHECK_TYPE}"
    )
  set(
    RUN_DIRECTORY
    "${EXECUTABLE_DIR_NAME}.${INPUT_FILE_NAME}.${CHECK_TYPE}"
    )
  if("${CHECK_TYPE}" STREQUAL "execute_check_output")
    set(_CHECK_OUTPUT_FILES "true")
  else()
    set(_CHECK_OUTPUT_FILES "false")
  endif()

  if ("${CHECK_TYPE}" STREQUAL "parse")
    add_test(
      NAME ${CTEST_NAME}
      COMMAND ${SPECTRE_TEST_RUNNER} ${CMAKE_BINARY_DIR}/bin/${EXECUTABLE}
      --check-options --input-file ${INPUT_FILE}
      )
  elseif("${CHECK_TYPE}" STREQUAL "execute" OR
         "${CHECK_TYPE}" STREQUAL "execute_check_output")
    add_test(
      NAME ${CTEST_NAME}
      # This script is written below, and only once
      COMMAND sh ${PROJECT_BINARY_DIR}/tmp/RunInputFileTest.sh
      ${EXECUTABLE} ${INPUT_FILE} ${RUN_DIRECTORY}
      ${EXPECTED_EXIT_CODE} ${_CHECK_OUTPUT_FILES}
      "${COMMAND_LINE_ARGS}"
      )
  else()
    message(FATAL_ERROR "Unknown check for input file: ${CHECK_TYPE}."
      "Known checks are: execute")
  endif()

  # Triple timeout if address sanitizer is enabled.
  if (ASAN)
    math(EXPR TIMEOUT "3 * ${TIMEOUT}")
  endif()

  spectre_test_timeout(TIMEOUT INPUT_FILE ${TIMEOUT})

  set_tests_properties(
    ${CTEST_NAME}
    PROPERTIES
    TIMEOUT ${TIMEOUT}
    LABELS "${TAGS}"
    ENVIRONMENT "${_INPUT_FILE_TEST_ENV_VARS}")
endfunction()

# Searches the directory INPUT_FILE_DIR for .yaml files and adds a test for each
# one. See `WritingTests.md` for details on controlling input file tests. Add
# input files to the whitelist at the bottom of this file to ignore those tests
function(add_input_file_tests INPUT_FILE_DIR INPUT_FILE_WHITELIST)
  set(INPUT_FILE_LIST "")
  file(GLOB_RECURSE INPUT_FILE_LIST ${INPUT_FILE_DIR} "${INPUT_FILE_DIR}*.yaml")
  set(TIMEOUT 2)
  list(TRANSFORM INPUT_FILE_WHITELIST PREPEND ${INPUT_FILE_DIR})

  foreach(INPUT_FILE ${INPUT_FILE_LIST})
    # Only parse the input file if we are allowed to
    if (${INPUT_FILE} IN_LIST INPUT_FILE_WHITELIST)
      continue()
    endif()
    if (INPUT_FILE MATCHES "\\.overlay_[0-9]*\\.yaml$")
      continue()
    endif()

    file(READ ${INPUT_FILE} INPUT_FILE_CONTENTS)

    # Read the priority of the test specified in input file, empty is accepted.
    string(REGEX MATCH "Priority:[^\n]+"
      INPUT_FILE_PRIORITY "${INPUT_FILE_CONTENTS}")
    if("${INPUT_FILE_PRIORITY}" STREQUAL "")
      set(INPUT_FILE_PRIORITY "normal")
    else()
      string(REGEX REPLACE "Priority:[ ]*" ""
        INPUT_FILE_PRIORITY "${INPUT_FILE_PRIORITY}")
      string(STRIP "${INPUT_FILE_PRIORITY}" INPUT_FILE_PRIORITY)
    endif()
    # Translate strings to numbers
    string(TOLOWER "${INPUT_FILE_PRIORITY}" INPUT_FILE_PRIORITY)
    if (${INPUT_FILE_PRIORITY} STREQUAL "high")
      set(INPUT_FILE_PRIORITY 1)
    elseif (${INPUT_FILE_PRIORITY} STREQUAL "normal")
      set(INPUT_FILE_PRIORITY 0)
    elseif (${INPUT_FILE_PRIORITY} STREQUAL "low")
      set(INPUT_FILE_PRIORITY -1)
    else()
      message(FATAL_ERROR "Unknown priority in input file ${INPUT_FILE}: "
        "${INPUT_FILE_PRIORITY}. Possible values are: high, normal, low")
    endif()

    # Only add tests with at least the priority specified in
    # `SPECTRE_INPUT_FILE_TEST_MIN_PRIORITY`
    if (${INPUT_FILE_PRIORITY} LESS
        ${SPECTRE_INPUT_FILE_TEST_MIN_PRIORITY})
      continue()
    endif()

    # Check if the executable name is present
    string(REGEX MATCH "Executable:[^\n]+"
      INPUT_FILE_EXECUTABLE "${INPUT_FILE_CONTENTS}")
    if("${INPUT_FILE_EXECUTABLE}" STREQUAL "")
      message(FATAL_ERROR "Could not find the executable in input "
        "file ${INPUT_FILE}. You must supply a line of the form:"
        "'Executable: EXECUTABLE_NAME'")
    endif()
    # Extract executable name, and remove trailing white space
    string(REGEX REPLACE "Executable:[ ]*" ""
      INPUT_FILE_EXECUTABLE "${INPUT_FILE_EXECUTABLE}")
    string(STRIP "${INPUT_FILE_EXECUTABLE}" INPUT_FILE_EXECUTABLE)

    string(REGEX MATCH "CommandLineArgs:[^\n]+"
      COMMAND_LINE_ARGS "${INPUT_FILE_CONTENTS}")
    string(REGEX REPLACE "CommandLineArgs:[ ]*" ""
      COMMAND_LINE_ARGS "${COMMAND_LINE_ARGS}")
    string(STRIP "${COMMAND_LINE_ARGS}" COMMAND_LINE_ARGS)

    string(REGEX MATCH "ExpectedExitCode:[^\n]+"
      EXPECTED_EXIT_CODE "${INPUT_FILE_CONTENTS}")
    string(REGEX REPLACE "ExpectedExitCode:[ ]*" ""
      EXPECTED_EXIT_CODE "${EXPECTED_EXIT_CODE}")
    string(STRIP "${EXPECTED_EXIT_CODE}" EXPECTED_EXIT_CODE)
    if("${EXPECTED_EXIT_CODE}" STREQUAL "")
      set(EXPECTED_EXIT_CODE "0")
    endif()

    # Read what tests to do. Currently "execute" and "parse" are available.
    string(REGEX MATCH "Check:[^\n]+"
      INPUT_FILE_CHECKS "${INPUT_FILE_CONTENTS}")
    # Extract list of checks to perform
    string(REGEX REPLACE "Check:[ ]*" ""
      INPUT_FILE_CHECKS "${INPUT_FILE_CHECKS}")
    string(STRIP "${INPUT_FILE_CHECKS}" INPUT_FILE_CHECKS)
    set(INPUT_FILE_CHECKS "${INPUT_FILE_CHECKS}")
    list(REMOVE_DUPLICATES "INPUT_FILE_CHECKS")
    # Convert all the checks to lower case to make life easier.
    string(TOLOWER "${INPUT_FILE_CHECKS}" INPUT_FILE_CHECKS)

    # Make sure that the 'parse' check is listed. If not, print an
    # error message that explains that it's needed and why.
    list(FIND "INPUT_FILE_CHECKS" "parse" FOUND_PARSE)
    if (${FOUND_PARSE} EQUAL -1)
      message(FATAL_ERROR
        "The input file: "
        "'${INPUT_FILE}' "
        "does not specify the 'parse' check. All input file tests must"
        " specify the 'parse' check which runs the executable passing"
        " the '--check-options' flag. With this flag the executable"
        " should check that the input file parses correctly and that"
        " the values specified in the input file do not violate any"
        " bounds or sanity checks.")
    endif (${FOUND_PARSE} EQUAL -1)

    # Read the timeout duration specified in input file, empty is accepted.
    # The default duration is 2 seconds.
    string(REGEX MATCH "Timeout:[^\n]+"
      INPUT_FILE_TIMEOUT "${INPUT_FILE_CONTENTS}")
    if("${INPUT_FILE_TIMEOUT}" STREQUAL "")
      set(INPUT_FILE_TIMEOUT "${TIMEOUT}")
    else()
      string(REGEX REPLACE "Timeout:[ ]*" ""
        INPUT_FILE_TIMEOUT "${INPUT_FILE_TIMEOUT}")
      string(STRIP "${INPUT_FILE_TIMEOUT}" INPUT_FILE_TIMEOUT)
    endif()

    foreach(CHECK_TYPE ${INPUT_FILE_CHECKS})
      add_single_input_file_test(
        ${INPUT_FILE}
        ${INPUT_FILE_EXECUTABLE}
        "${COMMAND_LINE_ARGS}"
        ${CHECK_TYPE}
        ${INPUT_FILE_TIMEOUT}
        ${EXPECTED_EXIT_CODE}
        )
    endforeach()
    add_dependencies(test-executables ${INPUT_FILE_EXECUTABLE})
  endforeach()
endfunction()

# Dependencies will be added as the tests are processed.
add_custom_target(test-executables)

# Write command to execute an input file and clean its output into a shell
# script, which makes it easier to chain multiple commands
configure_file(
  ${CMAKE_SOURCE_DIR}/cmake/RunInputFileTest.sh
  ${PROJECT_BINARY_DIR}/tmp/RunInputFileTest.sh
  @ONLY)

# These paths should be relative to the input file directory passed to
# `add_input_file_tests`
set(INPUT_FILE_WHITELIST
    "PreprocessCceWorldtube/PreprocessCceWorldtube.yaml")

add_input_file_tests("${CMAKE_SOURCE_DIR}/tests/InputFiles/" ${INPUT_FILE_WHITELIST})
