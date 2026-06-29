find_package(MPI COMPONENTS CXX)

if (MPI_CXX_FOUND)
    message(STATUS "Found MPI: ${MPI_CXX_COMPILER}")
    target_link_libraries(fptools PRIVATE MPI::MPI_CXX)
    target_compile_definitions(fptools PRIVATE USE_MPI)
else ()
    message(FATAL_ERROR
        "FPTOOLS_USE_MPI is ON but MPI was not found.\n"
        "   Linux: apt install libopenmpi-dev\n"
        "   macOS: brew install open-mpi"
        "   HPC:   module load openmpi (or the appropriate MPI module)\n"
    )
endif ()
