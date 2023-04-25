cmake_minimum_required(VERSION 3.10 FATAL_ERROR)

project(AS3 VERSION 1.0)

add_definitions(-std=c++20)
set(CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS, "${CXX_FLAGS}")
set(CMAKE_BUILD_TYPE Debug)
find_package(implot CONFIG REQUIRED)

#include_directories(matplotplusplus/include)
add_executable(main main.cpp)
target_link_libraries(main PRIVATE implot::implot)
#add_subdirectory(matplotplusplus)
#target_link_libraries(main PUBLIC matplot)
#find_package(glfw3 3.3 REQUIRED)
#find_package(OpenGL REQUIRED)

#target_include_directories(main PUBLIC ${OPENGL_INCLUDE_DIR})

#target_link_libraries(main "-framework Cocoa")
#target_link_libraries(main "-framework OpenGL")
#target_link_libraries(main "-framework IOKit")
#target_link_libraries(main glfw ${OPENGL_gl_LIBRARY})