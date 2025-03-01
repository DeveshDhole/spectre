# Distributed under the MIT License.
# See LICENSE.txt for details.

# Have to compile a file (can be empty) to check for flags.
# Kokkos' nvcc_wrapper doesn't support taking the file as stdin, so we
# have to write it to a file.
set(_CHECK_CXX_FLAGS_SOURCE "${CMAKE_BINARY_DIR}/CMakeFiles/CheckCxxFlags.cpp")
write_file(${_CHECK_CXX_FLAGS_SOURCE} "")

# Checks if a flag is supported by the compiler and creates the target
# TARGET_NAME whose INTERFACE_COMPILE_OPTIONS are set to the FLAG_TO_CHECK
# - LANGUAGE: language to check, setting the compiler and generated property
# - XTYPE: language as passed to the -x compiler flag
# - FLAG_TO_CHECK: the CXX flag to add if the compiler supports it
# - TARGET_NAME: the name of the target whose INTERFACE_COMPILE_OPTIONS are
#                set
function(create_compile_flag_target LANGUAGE XTYPE FLAG_TO_CHECK TARGET_NAME)
  # In order to check for a -Wno-* flag in gcc, you have to check the
  # -W* version instead.  See http://gcc.gnu.org/wiki/FAQ#wnowarning
  string(REGEX REPLACE ^-Wno- -W POSITIVE_FLAG_TO_CHECK ${FLAG_TO_CHECK})
  # Escape quotes for compiler command
  string(REPLACE "\"" "\\\"" POSITIVE_FLAG_TO_CHECK ${POSITIVE_FLAG_TO_CHECK})
  execute_process(
    COMMAND
    bash -c
    "LC_ALL=POSIX ${CMAKE_${LANGUAGE}_COMPILER} -Werror \
${POSITIVE_FLAG_TO_CHECK} -x ${XTYPE} \
-c ${_CHECK_CXX_FLAGS_SOURCE} -o /dev/null"
    RESULT_VARIABLE RESULT
    ERROR_VARIABLE ERROR_FROM_COMPILATION
    OUTPUT_QUIET)
  if(NOT TARGET ${TARGET_NAME})
    add_library(${TARGET_NAME} INTERFACE)
  endif(NOT TARGET ${TARGET_NAME})
  if(${RESULT} EQUAL 0)
    string(REPLACE " " ";" FLAG_TO_CHECK ${FLAG_TO_CHECK})
    set_property(TARGET ${TARGET_NAME}
      APPEND PROPERTY
      INTERFACE_COMPILE_OPTIONS
      $<$<COMPILE_LANGUAGE:${LANGUAGE}>:${FLAG_TO_CHECK}>)
  endif(${RESULT} EQUAL 0)
endfunction()

# Checks if a CXX flag is supported by the compiler and creates the target
# TARGET_NAME whose INTERFACE_COMPILE_OPTIONS are set to the FLAG_TO_CHECK
# - FLAG_TO_CHECK: the CXX flag to add if the compiler supports it
# - TARGET_NAME: the name of the target whose INTERFACE_COMPILE_OPTIONS are
#                set
function(create_cxx_flag_target FLAG_TO_CHECK TARGET_NAME)
  create_compile_flag_target(CXX c++ "${FLAG_TO_CHECK}" "${TARGET_NAME}")
endfunction()

# Same, but for C.
function(create_c_flag_target FLAG_TO_CHECK TARGET_NAME)
  create_compile_flag_target(C c "${FLAG_TO_CHECK}" "${TARGET_NAME}")
endfunction()

# Checks which of the CXX FLAGS_TO_CHECK are supported by the compiler
# and creates the target TARGET_NAME whose INTERFACE_COMPILE_OPTIONS
# are set to the FLAGS_TO_CHECK that are supported. If adding many flags,
# this will be much faster than calling create_cxx_flags_target multiple times.
# - LANGUAGE: language to check, setting the compiler and generated property
# - XTYPE: language as passed to the -x compiler flag
# - FLAGS_TO_CHECK: a semicolon separated string of CXX flags to try to add
#                   for compilation.
# - TARGET_NAME: the name of the target whose INTERFACE_COMPILE_OPTIONS are
#                set
function(create_compile_flags_target LANGUAGE XTYPE FLAGS_TO_CHECK TARGET_NAME)
  # In order to check for a -Wno-* flag in gcc, you have to check the
  # -W* version instead.  See http://gcc.gnu.org/wiki/FAQ#wnowarning
  set(POSITIVE_FLAGS_TO_CHECK)
  foreach(FLAG_TO_CHECK ${FLAGS_TO_CHECK})
    string(REGEX REPLACE ^-Wno- -W POSITIVE_FLAG_TO_CHECK ${FLAG_TO_CHECK})
    list(APPEND POSITIVE_FLAGS_TO_CHECK ${POSITIVE_FLAG_TO_CHECK})
  endforeach()
  string(REPLACE ";" " "
    POSITIVE_FLAGS_WITH_SPACES "${POSITIVE_FLAGS_TO_CHECK}")
  execute_process(
    COMMAND
    bash -c
    "LC_ALL=POSIX ${CMAKE_${LANGUAGE}_COMPILER} -Werror \
${POSITIVE_FLAGS_WITH_SPACES} -x ${XTYPE} \
-c ${_CHECK_CXX_FLAGS_SOURCE} -o /dev/null"
    RESULT_VARIABLE RESULT
    ERROR_VARIABLE ERROR_FROM_COMPILATION
    OUTPUT_QUIET)

  if(NOT TARGET ${TARGET_NAME})
    add_library(${TARGET_NAME} INTERFACE)
  endif(NOT TARGET ${TARGET_NAME})
  if(${RESULT} EQUAL 0)
    set_property(TARGET ${TARGET_NAME}
      APPEND PROPERTY
      INTERFACE_COMPILE_OPTIONS
      $<$<COMPILE_LANGUAGE:${LANGUAGE}>:${FLAGS_TO_CHECK}>)
  else(${RESULT} EQUAL 0)
    # Check each flag to see if it was marked as "invalid" in the output
    unset(FLAGS_TO_ADD)
    foreach(FLAG ${POSITIVE_FLAGS_TO_CHECK})
      string(FIND "${ERROR_FROM_COMPILATION}" "'${FLAG}'" FOUND_POS)
      if(${FOUND_POS} EQUAL -1)
        # For some reason:
        # list(FIND ${POSITIVE_FLAGS_TO_CHECK} ${FLAG} INDEX_OF_FLAG)
        # doesn't work with some compilers. This makes no sense but such is
        # life. As a work around we basically implement a find manually.

        # Find the index of the current flag in the POSITIVE_FLAGS_TO_CHECK
        # list. This is the index we use to get the original flag in the
        # FLAGS_TO_CHECK list.
        set(INDEX 0)
        foreach(POS_FLAG ${POSITIVE_FLAGS_TO_CHECK})
          if("${POS_FLAG}" STREQUAL "${FLAG}")
            break()
          endif()
          MATH(EXPR INDEX "${INDEX}+1")
        endforeach()
        set(TARGET_INDEX ${INDEX})
        set(INDEX 0)
        # Get original flag
        set(NEW_FLAG "")
        foreach(ORIGINAL_FLAG ${FLAGS_TO_CHECK})
          if(${INDEX} EQUAL ${TARGET_INDEX})
            set(NEW_FLAG ${ORIGINAL_FLAG})
            break()
          endif()
          MATH(EXPR INDEX "${INDEX}+1")
        endforeach()
        # Add the flag to the list of flags to add.
        set(FLAGS_TO_ADD "${FLAGS_TO_ADD};${NEW_FLAG}")
      endif(${FOUND_POS} EQUAL -1)
    endforeach(FLAG ${FLAGS_TO_CHECK})
    set_property(TARGET ${TARGET_NAME}
      APPEND PROPERTY
      INTERFACE_COMPILE_OPTIONS
      $<$<COMPILE_LANGUAGE:${LANGUAGE}>:${FLAGS_TO_ADD}>)

  endif(${RESULT} EQUAL 0)
endfunction()

# Checks which of the CXX FLAGS_TO_CHECK are supported by the compiler
# and creates the target TARGET_NAME whose INTERFACE_COMPILE_OPTIONS
# are set to the FLAGS_TO_CHECK that are supported. If adding many flags,
# this will be much faster than calling create_cxx_flags_target multiple times.
# - FLAGS_TO_CHECK: a semicolon separated string of CXX flags to try to add
#                   for compilation.
# - TARGET_NAME: the name of the target whose INTERFACE_COMPILE_OPTIONS are
#                set
function(create_cxx_flags_target FLAGS_TO_CHECK TARGET_NAME)
  create_compile_flags_target(CXX c++ "${FLAGS_TO_CHECK}" "${TARGET_NAME}")
endfunction()

# Same, but for C.
function(create_c_flags_target FLAGS_TO_CHECK TARGET_NAME)
  create_compile_flags_target(C c "${FLAGS_TO_CHECK}" "${TARGET_NAME}")
endfunction()

set(CMAKE_SUPPORTS_LINK_OPTIONS OFF)
if(CMAKE_VERSION VERSION_EQUAL 3.13 OR CMAKE_VERSION VERSION_GREATER 3.13)
  set(CMAKE_SUPPORTS_LINK_OPTIONS ON)
endif(CMAKE_VERSION VERSION_EQUAL 3.13 OR CMAKE_VERSION VERSION_GREATER 3.13)

if(CMAKE_SUPPORTS_LINK_OPTIONS)
  # Creates a target named TARGET_NAME that, if the linker flag FLAG_TO_CHECK
  # is supported, defines ${FLAG_TO_CHECK} as an INTERFACE_LINK_OPTION
  function(create_cxx_link_flag_target FLAG_TO_CHECK TARGET_NAME)
    include(CheckLinkerFlag)
    unset(CXX_LINKER_FLAG_WORKS CACHE)
    set(CMAKE_REQUIRED_QUIET 1)
    check_linker_flag(CXX ${FLAG_TO_CHECK} CXX_LINKER_FLAG_WORKS)
    unset(CMAKE_REQUIRED_QUIET)

    add_library(${TARGET_NAME} INTERFACE)
    if(CXX_LINKER_FLAG_WORKS)
      set_property(TARGET ${TARGET_NAME}
        APPEND PROPERTY
        INTERFACE_LINK_OPTIONS $<$<COMPILE_LANGUAGE:CXX>:${FLAG_TO_CHECK}>)
    endif()
  endfunction()
endif(CMAKE_SUPPORTS_LINK_OPTIONS)

# Checks if a flag is supported by the linker and adds it if it is
function(check_and_add_cxx_link_flag FLAG_TO_CHECK)
  include(CheckLinkerFlag)
  unset(CXX_LINKER_FLAG_WORKS CACHE)
  set(CMAKE_REQUIRED_QUIET 1)
  check_linker_flag(CXX ${FLAG_TO_CHECK} CXX_LINKER_FLAG_WORKS)
  unset(CMAKE_REQUIRED_QUIET)
  if(CXX_LINKER_FLAG_WORKS)
    set(CMAKE_CXX_LINK_FLAGS
      "${CMAKE_CXX_LINK_FLAGS} ${FLAG_TO_CHECK}" PARENT_SCOPE)
  endif()
endfunction()
