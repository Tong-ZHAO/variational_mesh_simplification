cmake_minimum_required(VERSION 3.9)

project(MeshSimplify)

# Instruct CMake to run moc automatically when needed.
set(CMAKE_AUTOMOC ON)

# Test if compiler has C++17
list(FIND CMAKE_CXX_COMPILE_FEATURES cxx_generalized_initializers has_cpp17)
if(has_cpp17 LESS 0)
  message(STATUS "NOTICE: These tests require a C++17 compiler and will not be compiled.")
  return()
endif()

if(COMMAND cmake_policy )
  cmake_policy(SET CMP0003 NEW)  
endif()

set(CMAKE_CXX_STANDARD 17)

find_package(CGAL REQUIRED)
find_package(Eigen3 3.1.0)

file(GLOB_RECURSE INCS "*.h")

if(CGAL_FOUND AND EIGEN3_FOUND)
    include_directories(BEFORE . ./include)
    include(${CGAL_USE_FILE})

    SET(MY_BUILD_PATH ${CMAKE_CURRENT_SOURCE_DIR}/build) #Variable build path
    SET(MY_SRC_PATH ${CMAKE_CURRENT_SOURCE_DIR})  #Variable src path
    INCLUDE_DIRECTORIES(${MY_BUILD_PATH})
    INCLUDE_DIRECTORIES(${MY_SRC_PATH})

    add_executable(MeshSimplify main.cpp ${INCS})

    target_link_libraries(MeshSimplify CGAL::CGAL)
else()
  message(STATUS "NOTICE: This program requires CGAL and Eigen. It will not be compiled.")

endif()