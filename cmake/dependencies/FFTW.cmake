find_package(PkgConfig QUIET)

# PRIMARY PATH: use pkg-config for clean IMPORTED target support.
if (PkgConfig_FOUND)
    pkg_check_modules(FFTW3 IMPORTED_TARGET fftw3)
    if (FFTW3_FOUND)
        message(STATUS "Found FFTW3 via pkg-config: ${FFTW3_PREFIX}")
        target_link_libraries(fptools PRIVATE PkgConfig::FFTW3)
        return()
    endif ()
endif ()

# FALLBACK PATH: perform a manual search.
find_library(FFTW3_LIB
    NAMES fftw3
    HINTS
        ${FFTW_ROOT}/lib
        $ENV{FFTW_ROOT}/lib
        $ENV{FFTW_DIR}/lib
        /usr/local/lib
        /opt/homebrew/lib
        /usr/lib
)

find_path(FFTW3_INCLUDE_DIR
    NAMES fftw3.h
    HINTS
        ${FFTW_ROOT}/include
        $ENV{FFTW_ROOT}/include
        $ENV{FFTW_DIR}/include
        /usr/local/include
        /opt/homebrew/include
        /usr/include
)

if (FFTW3_LIB AND FFTW3_INCLUDE_DIR)
    message(STATUS "Found FFTW3 (fallback): ${FFTW3_LIB}")
    target_include_directories(fptools PRIVATE ${FFTW3_INCLUDE_DIR})
    target_link_libraries(fptools PRIVATE ${FFTW3_LIB})
else ()
    message(FATAL_ERROR
        "FFTW3 not found. Please install FFTW3 and try again.\n"
        "   Linux: apt install libfftw3-dev\n"
        "   macOS: brew install fftw\n"
        "   HPC:   module load fftw\n"
        "If FFTW is installed in a non-standard location, set:\n"
        "   -DFFTW_ROOT=/path/to/fftw"
    )
endif ()
