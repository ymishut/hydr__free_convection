cmake_minimum_required(VERSION 2.8)
add_compile_options(-std=c++11 -Wall -pthread)

set(PROJECT_NAME free_conv)
project(${PROJECT_NAME})

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)
find_package(Boost REQUIRED)
find_package(Armadillo REQUIRED)

set(
  hs_sources
  src/animation.cpp
  src/basicflow.cpp
  src/calculation.cpp
  src/filereading.cpp
  src/galerk_method.cpp
  src/initialization.cpp
  src/inputdata_by_file.cpp
  src/laguerre_polynomials.cpp
  src/maindatastorage.cpp
  src/nachts_method.cpp
  src/stabilitypoints.cpp
  src/wide_output.cpp
)

add_executable(
  ${PROJECT_NAME}
  ${hs_sources}
  src/main.cpp
)

include_directories(${OPENGL_INCLUDE_DIRS} ${GLUT_INCLUDE_DIRS} ${BOOST_INCLUDE_DIRS} ${ARMADILLO_INCLUDE_DIRS})

target_link_libraries(${PROJECT_NAME} ${OPENGL_LIBRARIES} ${GLUT_LIBRARY} ${ARMADILLO_LIBRARIES} pthread)
