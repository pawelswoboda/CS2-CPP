cmake_minimum_required(VERSION 2.8.12)

# C++11
add_compile_options(-std=c++14)

#SET(C++_STD_FLAG "c++14")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=${C++_STD_FLAG}")
#MARK_AS_ADVANCED(C++_STD_FLAG)

# compiler options
add_definitions(-DIL_STD)
add_definitions(-ffast-math)
add_definitions(-march=native)
#add_definitions(-libc++)
add_definitions(-libstdc++)

option(BUILD_TESTS "build tests" OFF)

if(BUILD_TESTS)
   option(DOWNLOAD_TEST_DEPENDENCIES "Download test dependencies" ON)
   # download external projects here
   if(DOWNLOAD_TEST_DEPENDENCIES)
      project(DOWNLOAD_TEST_DEPENDENCIES_PROJECT NONE)

      include(ExternalProject)
      set_property(DIRECTORY PROPERTY EP_BASE Dependencies)
      set(DEPENDENCIES)
      set(EXTRA_CMAKE_ARGS)

      # unit tests framework
      list(APPEND DEPENDENCIES Catch_Project)
      ExternalProject_Add(
         Catch_Project
         GIT_REPOSITORY "https://github.com/philsquared/Catch.git"
         GIT_TAG "v1.5.6"
         UPDATE_COMMAND ""
         INSTALL_COMMAND ""
         BUILD_COMMAND ""
         CONFIGURE_COMMAND ""
         )
      ExternalProject_Get_Property(Catch_Project source_dir)
      set(CATCH_INCLUDE_DIR ${source_dir}/include CACHE INTERNAL "Path to include folder for Catch")

      ExternalProject_Add(CS2_CPP
         DEPENDS ${DEPENDENCIES}
         SOURCE_DIR ${PROJECT_SOURCE_DIR}
         CMAKE_ARGS -DDOWNLOAD_TEST_DEPENDENCIES=OFF ${EXTRA_CMAKE_ARGS}
         INSTALL_COMMAND ""
         BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR})

      return()
   endif()
endif()

project(CS2_CPP)
#add_library(CS2_CPP STATIC
#   mcmf.cpp mcmf.h)

add_subdirectory(test)

