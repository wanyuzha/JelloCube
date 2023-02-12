cmake_minimum_required(VERSION 3.8)
set (CMAKE_CXX_STANDARD 11)

project("hw1")

find_package(OpenGL REQUIRED)
find_package(GLUT REQUIRED)

include_directories(
        ${OPENGL_INCLUDE_DIRS}
        ${GLUT_INCLUDE_DIRS}
)

add_executable(jello
        "jello.h"
        "jello.cpp"
        "input.h"
        "input.cpp"
        "openGL-headers.h"
        "physics.h"
        "physics.cpp"
        "pic.h"
        "pic.cpp"
        "ppm.cpp"
        "showCube.h"
        "showCube.cpp"
        )

target_link_libraries(jello
        OpenGL::GL
        GLUT::GLUT
        )
