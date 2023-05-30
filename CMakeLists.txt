cmake_minimum_required(VERSION 3.1)
project(deformation)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

### libIGL options: choose between header only and compiled static library
option(LIBIGL_USE_STATIC_LIBRARY     "Use libigl as static library" ON)
option(LIBIGL_GLFW       "Use GLFW"           ON)

include(libigl)

# Add your project files
include_directories("include/")
include_directories("E:/dev/libigl-2.4.0/build/install/include/")
include_directories("E:/dev/libigl-2.4.0/include/igl/")
include_directories("E:/dev/libigl-2.4.0/include/")
include_directories("E:/dev/glad_3.3_core/include/")
include_directories("E:/dev/glfw-3.3.7/build/install/include/")
include_directories("E:/dev/eigen-3.3.9/")
link_directories("E:/dev/libigl-2.4.0/build/install/lib_release/")
link_directories("E:/dev/glfw-3.3.7/build/install/lib/")

if(USE_SOLUTION)
  file(GLOB SRCFILES solution/*.cpp)
else()
  file(GLOB SRCFILES source/*.cpp)
endif()

add_library(core ${SRCFILES})

add_executable(deformation "main.cpp")
target_link_libraries(deformation core opengl32 glfw3)
