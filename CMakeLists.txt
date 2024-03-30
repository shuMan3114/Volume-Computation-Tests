project( VolEsti-cpp-example )

CMAKE_MINIMUM_REQUIRED(VERSION 3.11)

set(CMAKE_ALLOW_LOOSE_LOOP_CONSTRUCTS true)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

add_definitions(-DDISABLE_NLP_ORACLES)

option(BUILTIN_EIGEN "Use eigen from /volesti/external" OFF)

include("/volesti/external/cmake-files/Eigen.cmake")
GetEigen()

include("/volesti/external/cmake-files/Boost.cmake")
GetBoost()

include("/volesti/external/cmake-files/LPSolve.cmake")
GetLPSolve()

# Find lpsolve library
find_library(LP_SOLVE NAMES liblpsolve55.so PATHS /usr/lib/lp_solve)

if (NOT LP_SOLVE)
  message(FATAL_ERROR "This program requires the lp_solve library, and will not be compiled.")
else ()
  message(STATUS "Library lp_solve found: ${LP_SOLVE}")

  set(CMAKE_EXPORT_COMPILE_COMMANDS "ON")

  include_directories (BEFORE /volesti/external)
  include_directories (BEFORE /volesti/include)
  include_directories (BEFORE /doctest/doctest)

  # for Eigen
  if (${CMAKE_VERSION} VERSION_LESS "3.12.0")
    add_compile_options(-D "EIGEN_NO_DEBUG")
  else ()
    add_compile_definitions("EIGEN_NO_DEBUG")
  endif ()


  add_definitions(${CMAKE_CXX_FLAGS} "-std=c++11")  # enable C++11 standard

  add_executable (volume_estimation_testing_vpoly volume_cg_vpolytope.cpp)
  TARGET_LINK_LIBRARIES(volume_estimation_testing_vpoly ${LP_SOLVE})

endif()
