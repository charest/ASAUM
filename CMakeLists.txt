#------------------------------------------------------------------------------#
# Set the minimum CMake version
#------------------------------------------------------------------------------#

cmake_minimum_required(VERSION 3.18)

#------------------------------------------------------------------------------#
# Setup the project
#------------------------------------------------------------------------------#

project(ASAUM VERSION 0.1.0)

set(default_build_type "Debug")
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message(STATUS "Setting build type to '${default_build_type}' as none was specified.")
  set(CMAKE_BUILD_TYPE "${default_build_type}" CACHE
      STRING "Choose the type of build." FORCE)
  # Set the possible values of build type for cmake-gui
  set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS
    "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()


# cmake module path
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake")

# We need C++ 17
add_library( asaum )
set_property(TARGET asaum PROPERTY CXX_STANDARD 14)
set_property(TARGET asaum PROPERTY CXX_STANDARD_REQUIRED on)
set_property(TARGET asaum PROPERTY CXX_EXTENSIONS off)

include(CTest) # note: this adds a BUILD_TESTING which defaults to ON

set( ASAUM_TEST_TOLERANCE 2e-14 CACHE INTERNAL "")
set( ASAUM_TOOL_DIR "${PROJECT_SOURCE_DIR}/tools" CACHE INTERNAL "")

#------------------------------------------------------------------------------#
# Setup third party packages
#------------------------------------------------------------------------------#

find_package(Threads QUIET)
option(HAVE_THREADS "Have the pthreads library" ${Threads_FOUND})
if (HAVE_THREADS)
  target_link_libraries(asaum PUBLIC Threads::Threads)
endif()

find_package(EXODUSII QUIET)
option(HAVE_EXODUS "Have the exodus library" ${EXODUSII_FOUND})
if(HAVE_EXODUS)
  target_include_directories(asaum PUBLIC ${EXODUSII_INCLUDE_DIRS})
  target_link_libraries(asaum PUBLIC ${EXODUSII_LIBRARIES})
else()
#  set_target_properties(exogen PROPERTIES EXCLUDE_FROM_ALL true)
endif()

find_package(MPI COMPONENTS C CXX)
if (MPI_FOUND)
  target_link_libraries(asaum PUBLIC MPI::MPI_C MPI::MPI_CXX)
endif()


find_package(ParMetis 4.0 QUIET)
option(HAVE_PARMETIS "Have the parmetis library" ${ParMetis_FOUND})
if (HAVE_PARMETIS)

  find_package(METIS 5.1)
  target_include_directories(asaum PUBLIC ${METIS_INCLUDE_DIRS})
  target_link_libraries(asaum PUBLIC ${METIS_LIBRARIES})

  target_include_directories(asaum PUBLIC ${ParMetis_INCLUDE_DIRS})
  target_link_libraries(asaum PUBLIC ${ParMetis_LIBRARIES})
endif()

find_package(Lua 5.4 QUIET)
if (NOT Lua_FOUND)
  find_package(Lua 5.2 REQUIRED)
endif()
target_include_directories(asaum PUBLIC ${LUA_INCLUDE_DIR})
target_link_libraries(asaum PUBLIC ${LUA_LIBRARIES} ${CMAKE_DL_LIBS})


# vtk
find_package(ZLIB QUIET)
option(HAVE_ZLIB "Have the zlib library" ${ZLIB_FOUND})
if (HAVE_ZLIB)
  target_link_libraries(asaum PUBLIC ZLIB::ZLIB)
endif()


option(ENABLE_FP_EXCEPTIONS "Enable GNU floating point exceptions" off)
option(ENABLE_TRACE "Enable tracing" off)



#------------------------------------------------------------------------------#
# configure library
#------------------------------------------------------------------------------#

add_subdirectory(libsrc)
add_subdirectory(example)
#if (HAVE_EXODUS)
#  add_subdirectory(exogen)
#endif()

#------------------------------------------------------------------------------#
# Add testing
#------------------------------------------------------------------------------#
if (BUILD_TESTING)
  include(regression)

  # find python for running regression tests
  find_package (Python COMPONENTS Interpreter REQUIRED)
  
  add_subdirectory(tests)
endif()

