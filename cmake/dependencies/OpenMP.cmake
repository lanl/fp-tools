find_package(OpenMP)

if (OpenMP_CXX_FOUND)
    message(STATUS "Found OpenMP ${OpenMP_CXX_VERSION} -- enabling OpenMP support")
    target_link_libraries(fptools PRIVATE OpenMP::OpenMP_CXX)
    target_compile_definitions(fptools PRIVATE USE_OMP)
else ()
    message(WARNING
        "OpenMP is required but was not found.\n"
        "   Linux: install libomp-dev "
        "   macOS: brew install libomp\n"
    )
endif ()
