cmake_minimum_required(VERSION 3.13)
project(pdc_mini_aevol)

set(CMAKE_CXX_STANDARD 14)
set(TARGET "micro_aevol")
list(APPEND CMAKE_MODULE_PATH "${KOKKOS_BUILD_DIR}")
list(APPEND CMAKE_PREFIX_PATH "${KOKKOS_INSTALL_DIR}")
list(APPEND CMAKE_PREFIX_PATH "${KOKKOS_INSTALL_DIR}/lib/CMake/Kokkos")
list(APPEND CMAKE_PREFIX_PATH "${KOKKOS_INSTALL_DIR}/lib/cmake/Kokkos")

set(DO_TRACES OFF CACHE BOOL "Whether to enable the Traces library")
# set(USE_OMP OFF CACHE BOOL "Whether to enable OpenMP parallelization")
# set(USE_CUDA OFF CACHE BOOL "Whether to enable CUDA parallelization")

if ( DO_TRACES )
    add_definitions(-DTRACES)
    message( STATUS "Traces are activated" )
endif ( DO_TRACES )

# if ( USE_CUDA )
#     # https://cliutils.gitlab.io/modern-cmake/chapters/packages/CUDA.html
#     find_package(CUDA REQUIRED)
#     include(CheckLanguage)
#     check_language(CUDA)
#     if(CMAKE_CUDA_COMPILER)
#         enable_language(CUDA)
#         set(CMAKE_CUDA_STANDARD 14)
#         set(CMAKE_CUDA_STANDARD_REQUIRED ON)
#         # https://arnon.dk/tag/nvcc-flags/
#         set(CMAKE_CUDA_FLAGS "-arch=sm_60")
#         set(CMAKE_CUDA_FLAGS_DEBUG "-Xcompiler -rdynamic -lineinfo")
#     else()
#         message(FATAL_ERROR "No CUDA support found")
#     endif()
# endif ( USE_CUDA )

# if ( USE_OMP )
#     find_package(OpenMP REQUIRED)
# endif ( USE_OMP )

find_package(ZLIB REQUIRED)
find_package(Kokkos REQUIRED)



# add_library(
#     Abstract_ExpManager.cpp
#     ExpManager.cpp
#     AeTime.cpp
#     DnaMutator.cpp
#     MutationEvent.cpp
#     Organism.cpp
#     Stats.cpp
#     Threefry.cpp
#     Dna.cpp
# )

add_executable(${TARGET} main.cpp)
target_include_directories(${TARGET} PUBLIC "${KOKKOS_INSTALL_DIR}/include")
target_link_libraries(${TARGET} PUBLIC ZLIB::ZLIB Kokkos::kokkos
    Abstract_ExpManager.cpp
    ExpManager.cpp
    AeTime.cpp
    DnaMutator.cpp
    MutationEvent.cpp
    Organism.cpp
    Stats.cpp
    Threefry.cpp
    Dna.cpp
)


# if ( OPENMP_FOUND )
#     target_link_libraries(micro_aevol PUBLIC OpenMP::OpenMP_CXX)
#     target_compile_definitions(micro_aevol PUBLIC USE_OMP)
# endif ()

# add_executable(${TARGET}_cpu main.cpp)
# target_link_libraries(micro_aevol_cpu ${TARGET})

# if ( USE_CUDA )
#     add_subdirectory(cuda)
#     add_executable(micro_aevol_gpu main.cpp)
#     # nvToolsExt for enhanced profiling (ad-hoc chunks)
#     target_link_libraries(micro_aevol_gpu PUBLIC cuda_micro_aevol micro_aevol)
# endif ()

# create a make clean rule for the project
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY ADDITIONAL_MAKE_CLEAN_FILES
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeCache.txt
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles
    ${CMAKE_CURRENT_BINARY_DIR}/cmake_install.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/Makefile
)