#===================================================================================================== PROJECT SETUP ===
cmake_minimum_required(VERSION 3.14)
project(tardigrade_micromorphic_tools)

# Set common project paths relative to project root directory
set(CPP_SRC_PATH "src/cpp")
set(CPP_TEST_PATH "${CPP_SRC_PATH}/tests")
set(CMAKE_SRC_PATH "src/cmake")
set(DOXYGEN_SRC_PATH "docs/doxygen")
set(SPHINX_SRC_PATH "docs/sphinx")
# Add the cmake folder to locate project CMake module(s)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/${CMAKE_SRC_PATH}" ${CMAKE_MODULE_PATH})

set(TARDIGRADE_MICROMORPHIC_TOOLS_BUILD_PYTHON_BINDINGS ON CACHE BOOL "Flag for whether the python bindings should be built")

# Add a flag for if a full build of all tardigrade repositories should be performed
set(TARDIGRADE_FULL_BUILD OFF CACHE BOOL "Flag for whether a full build of Tardigrade should be performed (i.e., all repos pulled from git and built)")

# Set build type checks
string(TOLOWER "${CMAKE_BUILD_TYPE}" cmake_build_type_lower)
set(upstream_required "")
set(not_conda_test "true")
set(project_link_string ${PROJECT_NAME})
if(cmake_build_type_lower STREQUAL "release")
    set(upstream_required "REQUIRED")
elseif(cmake_build_type_lower STREQUAL "conda-test")
    set(upstream_required "REQUIRED")
    set(not_conda_test "false")
    # Find the installed project library
    find_file(installed_linked_library
              "lib${PROJECT_NAME}${CMAKE_SHARED_LIBRARY_SUFFIX}"
              PATHS "$ENV{CONDA_PREFIX}/${CMAKE_INSTALL_LIBDIR}"
              NO_CACHE
              REQUIRED)
    set(project_link_string ${installed_linked_library})
endif()

# Get version number from Git
if(${not_conda_test} STREQUAL "true")
    if (DEFINED ENV{PREFIX})
        set(Python_ROOT_DIR "$ENV{PREFIX}/bin")
    endif()
    set(Python_FIND_STRATEGY LOCATION)
    find_package(Python COMPONENTS Interpreter REQUIRED)
    execute_process(COMMAND ${Python_EXECUTABLE} -m setuptools_scm
		            OUTPUT_VARIABLE ${PROJECT_NAME}_VERSION_STRING_FULL
                    OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(${${PROJECT_NAME}_VERSION_STRING_FULL} STREQUAL "")
        set(${PROJECT_NAME}_VERSION 0.0.0)
    else()
         string(REGEX MATCH "[0-9]+\.[0-9]+\.[0-9]+" ${PROJECT_NAME}_VERSION ${${PROJECT_NAME}_VERSION_STRING_FULL})
    endif()
    project(${PROJECT_NAME} VERSION ${${PROJECT_NAME}_VERSION})
endif()

# Add installation directory variables
include(GNUInstallDirs)

# Make the code position independent
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# Set the c++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic")
if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lrt -O3")
endif()

# Enable CTest
enable_testing()

#================================================================================================= FIND DEPENDENCIES ===
# Find eigen
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
include_directories(${EIGEN3_INCLUDE_DIR})
if(EIGEN3_FOUND)
    message(STATUS "Found Eigen3: ${EIGEN3_INCLUDE_DIR}")
endif()

# Find related, but external, projects in installed environment
include(FetchContent)
set(upstream_packages "tardigrade_error_tools" "tardigrade_vector_tools" "tardigrade_constitutive_tools")
if(${not_conda_test} STREQUAL "false")
    # During conda-build testing, we must find the installed project files as if they were an external project
    set(upstream_packages ${upstream_packages} ${PROJECT_NAME})
endif()
foreach(package ${upstream_packages})
    string(TOUPPER "${package}" package_upper)
    set(${package_upper}_BUILD_PYTHON_BINDINGS ${TARDIGRADE_MICROMORPHIC_TOOLS_BUILD_PYTHON_BINDINGS} CACHE INTERNAL "Setting ${package}'s python binding flag to the global value")
    if(TARDIGRADE_FULL_BUILD)
        message("Not attempting to find ${package}")
    else()
        find_package(${package} ${upstream_required} CONFIG)
    endif()
    if(${package}_FOUND)
        message(STATUS "Found ${package}: ${${package}_DIR}")
    else()
        # Find related, but external, projects using FetchContent and building locally
        message(WARNING "Did not find an installed ${package} package. Attempting local build with FetchContent.")
        if(NOT DEFINED ${package_upper}_FETCHCONTENT_VERSION)
            set(${package_upper}_FETCHCONTENT_VERSION "origin/dev")
        endif()
        message("${package_upper} is being built with version ${${package_upper}_FETCHCONTENT_VERSION}")
        FetchContent_Declare(
            ${package}
            GIT_REPOSITORY http://github.com/UCBoulder/${package}.git
            GIT_TAG        ${${package_upper}_FETCHCONTENT_VERSION}
        )
        FetchContent_MakeAvailable(${package})
    endif()
endforeach(package)

#=============================================================================================== ADD PROJECT TARGETS ===
# MUST COME AFTER DEPENDENCY LOCATING
# Add project source directories
if(${not_conda_test} STREQUAL "true")
    include_directories(${CPP_SRC_PATH})
    add_subdirectory(${CPP_SRC_PATH})
endif()

# Only add tests and documentation for current project builds. Protects downstream project builds.
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
    # Find Boost and add tests
    find_package(Boost 1.53.0 REQUIRED COMPONENTS unit_test_framework)
    add_subdirectory(${CPP_TEST_PATH})
    # Add docs
    if(${not_conda_test} STREQUAL "true")
        add_subdirectory(${DOXYGEN_SRC_PATH})
        add_subdirectory(${SPHINX_SRC_PATH})
    endif()
endif()

#==================================================================================== SETUP INSTALLATION CMAKE FILES ==
if(${not_conda_test} STREQUAL "true")
    foreach(package ${PROJECT_NAME})
        include(CMakePackageConfigHelpers)
        write_basic_package_version_file("${package}ConfigVersion.cmake"
                                         VERSION ${PROJECT_VERSION}
                                         COMPATIBILITY SameMajorVersion)
        configure_package_config_file(
          "${PROJECT_SOURCE_DIR}/${CMAKE_SRC_PATH}/Config.cmake.in"
          "${PROJECT_BINARY_DIR}/${package}Config.cmake"
          INSTALL_DESTINATION
          ${CMAKE_INSTALL_DATAROOTDIR}/${package}/cmake)

        # CMake won't build the targets for local builds of upstream projects
        if(cmake_build_type_lower STREQUAL release)
            install(EXPORT ${package}_Targets
                    FILE ${package}Targets.cmake
                    DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${package}/cmake)
        endif()

        install(FILES "${PROJECT_BINARY_DIR}/${package}Config.cmake"
                      "${PROJECT_BINARY_DIR}/${package}ConfigVersion.cmake"
                DESTINATION ${CMAKE_INSTALL_DATAROOTDIR}/${package}/cmake)
    endforeach(package)
endif()
