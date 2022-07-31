cmake_minimum_required(VERSION 3.10)

project(gcopter)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_BUILD_TYPE Release)

find_package(Python COMPONENTS Interpreter Development NumPy)
find_package(Eigen3 REQUIRED)

include(FetchContent)
FetchContent_Declare(yaml-cpp
  GIT_REPOSITORY https://github.com/jbeder/yaml-cpp.git)
FetchContent_GetProperties(yaml-cpp)

if(NOT yaml-cpp_POPULATED)
  message(STATUS "Fetching yaml-cpp...")
  FetchContent_Populate(yaml-cpp)
  add_subdirectory(${yaml-cpp_SOURCE_DIR} ${yaml-cpp_BINARY_DIR})
endif()

include_directories(
  include
  ${Python_INCLUDE_DIRS}
  ${Python_NumPy_INCLUDE_DIRS}
  ${EIGEN3_INCLUDE_DIRS}
)

add_executable(curve_gen src/curve_gen.cpp)
target_link_libraries(curve_gen PRIVATE ${Python_LIBRARIES} yaml-cpp::yaml-cpp)
