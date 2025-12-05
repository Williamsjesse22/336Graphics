cmake_minimum_required(VERSION 3.16)
project(raytracer)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_executable(raytracer
    src/main.cpp
    src/stb_image_impl.cpp
    src/tiny_obj_loader_impl.cpp
)

target_include_directories(raytracer PRIVATE src)
