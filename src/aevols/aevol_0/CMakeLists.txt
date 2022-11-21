cmake_minimum_required(VERSION 3.13)
project(pdc_mini_aevol)

set(CMAKE_CXX_STANDARD 14)
set(TARGET "micro_aevol")
set(TARGET_EXE ${TARGET}_kokkos)

list(APPEND CMAKE_MODULE_PATH $ENV{KOKKOS_BUILD_DIR})
list(APPEND CMAKE_PREFIX_PATH $ENV{KOKKOS_INSTALL_DIR}:)
list(APPEND CMAKE_PREFIX_PATH $ENV{KOKKOS_INSTALL_DIR}/lib/CMake/Kokkos:)
list(APPEND CMAKE_PREFIX_PATH $ENV{KOKKOS_INSTALL_DIR}/lib/cmake/Kokkos)
message("CMAKE_PREFIX_PATH = " ${CMAKE_PREFIX_PATH})

# preprocessor
set(DO_TRACES ON CACHE BOOL "Whether to enable the Traces library")

if ( DO_TRACES )
    add_definitions(-DTRACES)
    message( STATUS "Traces are activated" )
endif ( DO_TRACES )

# packages
find_package(ZLIB REQUIRED)
find_package(Kokkos REQUIRED)

# build TARGET as library
add_library(
    ${TARGET}
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
target_link_libraries(${TARGET} PUBLIC ZLIB::ZLIB Kokkos::kokkos)
target_include_directories(${TARGET} PUBLIC ${KOKKOS_INSTALL_DIR}/include)

# build TARGET_EXE as executable from TARGET library
add_executable(${TARGET_EXE} main.cpp)
target_link_libraries(${TARGET_EXE} ${TARGET})

# create a make clean rule for the project
set_property(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY ADDITIONAL_MAKE_CLEAN_FILES
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeCache.txt
    ${CMAKE_CURRENT_BINARY_DIR}/CMakeFiles
    ${CMAKE_CURRENT_BINARY_DIR}/cmake_install.cmake
    ${CMAKE_CURRENT_BINARY_DIR}/Makefile
)

# end
message("Done. " ${TARGET_EXE} " is ready to be executed.")
