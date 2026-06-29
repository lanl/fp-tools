
# fp_tools_set_default_build_type(<type>) #####################################
#
# Sets the default CMAKE_BUILD_TYPE if the user has not specified one. Also
# registers valid build types so CMake GUIs show a dropdown instead of a
# free-text field.
#
# Usage:
#   fp_tools_set_default_build_type("Release")
#
function(fp_tools_set_default_build_type default_build_type)
    set(_valid_build_types Release Debug RelWithDebInfo MinSizeRel)

    # Clean up empty strings or completely unset scenarios
    if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
        # 1. Force the global cache variable to change
        set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE STRING 
            "Build type (Release/Debug/RelWithDebInfo/MinSizeRel)" FORCE)
        
        # 2. Push it to the parent scope (for the rest of the CMake project)
        set(CMAKE_BUILD_TYPE "${default_build_type}" PARENT_SCOPE)
        
        # 3. Unset any local normal variable shadowing in this function scope
        unset(CMAKE_BUILD_TYPE) 
    endif ()

    # FALLBACK: if no local variable exists, CMake reads the CACHE variable.
    # To guarantee 'IN_LIST' evaluates the updated value, we pull it from the CACHE.
    if (NOT CMAKE_BUILD_TYPE)
        set(_current_build_type "$CACHE{CMAKE_BUILD_TYPE}")
    else()
        set(_current_build_type "${CMAKE_BUILD_TYPE}")
    endif()

    # Provide the selection strings for CMake GUI dropdowns
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS ${_valid_build_types})

    # Validate against the resolved value rather than the potentially shadowed variable
    if (NOT _current_build_type IN_LIST _valid_build_types)
        string(REPLACE ";" "; " _valid_build_types_printed "${_valid_build_types}")
        message(WARNING
            "CMAKE_BUILD_TYPE '${_current_build_type}' is not one of "
            "${_valid_build_types_printed}. Is this a typo?"
        )
    endif ()
endfunction ()


# fp_tools_set_cxx_warnings(<target>) #########################################
#
# Applies a standard set of compiler warnings to a target. Uses 
# target_compile_options so warnings are scoped to the target and do not
# leak into vendored dependencies.
#
# Covers GCC and Clang. MSVC is not a currently-supported platform for this 
# project, but a stub is included to avoid hard failures if someone tries.
#
# Usage:
#   fp_tools_set_cxx_warnings(fptools)
#
function (fp_tools_set_cxx_warnings target)
    if (CMAKE_CXX_COMPILER_ID MATCHES "GNU")
        target_compile_options(${target} PRIVATE
            -Wall                   # Enables a broad set of commonly useful warnings
            -Wextra                 # Enables additional warnings -Wall doesn't cover
            -Wpedantic              # Enforces strict ISO C++ compliance; catches non-standard extensions
            -Wshadow                # Warns when a variable declaration shadows an outer scope variable
            -Wunreachable-code      # Warms about code after a return statement or similar
            -Woverloaded-virtual    # Warns if a derived class accidentally hides a virtual method
            -Wno-comment            # Suppress multiline comment warnings in headers
        )
    elseif (CMAKE_CXX_COMPILER_ID MATCHES "Clang")
        target_compile_options(${target} PRIVATE
            -Wall
            -Wextra
            -Wpedantic
            -Wshadow
            -Wunreachable-code
            -Woverloaded-virtual
            -Wno-comment
        )
    elseif (CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
        # FP-Tools is not officially supported on Windows but we avoid hard
        # failures here in case someone tries.
        message(WARNING "MSVC is not an officially supported compiler for FP-Tools.")
        target_compile_options(${target} PRIVATE /W3)
    endif ()
endfunction ()

# fp_tools_print_summary() ####################################################
#
# Prints a human-readable build configuration summary at the end of CMake
# configuration. This is useful when building on unfamiliar HPC clusters
# to confirm that all dependencies were found correctlyand options are set as
# expected.
#
function (fp_tools_print_summary)
    message(STATUS "")
    message(STATUS "=============================================")
    message(STATUS "        FP-Tools Build Configuration")
    message(STATUS "=============================================")
    message(STATUS "    Version      : ${PROJECT_VERSION}")
    message(STATUS "    Build type   : ${CMAKE_BUILD_TYPE}")
    message(STATUS "    C++ standard : C++${CMAKE_CXX_STANDARD}")
    message(STATUS "    Compiler     : ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
    message(STATUS "    Install dir  : ${CMAKE_INSTALL_PREFIX}")
    message(STATUS "")
    message(STATUS "    Options:")
    message(STATUS "        MPI support     : ${FPTOOLS_USE_MPI}")
    message(STATUS "        OpenMP support  : ${OpenMP_CXX_FOUND}")
    message(STATUS "        Build tests     : ${FPTOOLS_BUILD_TESTS}")
    message(STATUS "")
    message(STATUS "    Output directories:")
    message(STATUS "        Binaries : ${CMAKE_RUNTIME_OUTPUT_DIRECTORY}")
    message(STATUS "        Libraries: ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
    message(STATUS "=============================================")
    message(STATUS "")
endfunction ()
