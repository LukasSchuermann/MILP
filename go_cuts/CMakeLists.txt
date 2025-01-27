cmake_minimum_required(VERSION 3.3)
set(CMAKE_CXX_STANDARD 17)

project(vrp)

if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif ()

set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -fopenmp -DNDEBUG")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")


### This Project #####
add_subdirectory(src)

set(SOURCES
        src/main_go.cpp
        src/ConshdlrGO.cpp
        src/event_setobjlimit.cpp
)


add_executable(vrp ${SOURCES})

target_include_directories(vrp
        PRIVATE
        ${PROJECT_SOURCE_DIR}/include
)

# for testing
#enable_testing()
#add_test(NAME vrpbuild/small-test.txt
#        COMMAND vrp ../instances/r101/r1_n20_nv25_p0.25_l10.json)
#add_test(vrpbuild/small-test.txt ${EXE})

# SCIP
# Note: SCIP_DIR as a path is only a recommendation for cmake
# if scip dir is not set
message(STATUS "Looking for SCIP.")
if (NOT SCIP_DIR)
    message(STATUS "No SCIP specified, will look for environment variable.")
    set(SCIP_DIR "$ENV{SCIP_DIR}" CACHE PATH "Path to SCIP directory")
else ()
    message(STATUS "Will use SCIP_DIR ${SCIP_DIR}.")
endif ()
find_package(SCIP PATHS ${SCIP_DIR} NO_DEFAULT_PATH REQUIRED)
message(STATUS "Using SCIP library from here ${SCIP_DIR}.")

## This project links

target_link_libraries(vrp
        ${SCIP_LIBRARIES}
)


# get_cmake_property(_variableNames VARIABLES)
# list (SORT _variableNames)
# foreach (_variableName ${_variableNames})
#     message(STATUS "${_variableName}=${${_variableName}}")
# endforeach()



