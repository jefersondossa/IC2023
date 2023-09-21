cmake_minimum_required (VERSION 3.10.0)

project (IC2023)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake")

# Enables NeoPZ library
find_package(NeoPZ REQUIRED HINTS ${CMAKE_SOURCE_DIR}/../install_neopz/ ${CMAKE_SOURCE_DIR}/install_neopz/ "../neopz_install")

message("Processor architecture = ${CMAKE_SYSTEM_PROCESSOR}")

# Project to apply error estimation for elasticity problems
add_executable(main main.cpp)
target_link_libraries(main PUBLIC NeoPZ::pz)

