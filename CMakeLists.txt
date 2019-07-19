cmake_minimum_required(VERSION 3.1)
SET(CMAKE_BUILD_TYPE Release)
project(specCoarsen)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  ON)

find_package(LIBIGL REQUIRED)
message("libigl location" ${LIBIGL_INCLUDE_DIR})

# Add your project files
include_directories("src/" "external/spectra/include/")

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "dir='${dir}'")
endforeach()

file(GLOB SRCFILES main.cpp src/*.cpp)
add_executable(${PROJECT_NAME}_bin ${SRCFILES})
target_link_libraries(${PROJECT_NAME}_bin igl::core igl::opengl igl::opengl_glfw)



