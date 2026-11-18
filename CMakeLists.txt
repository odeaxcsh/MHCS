cmake_minimum_required(VERSION 3.21)

project(PKMC LANGUAGES CXX)

add_library(PKMC STATIC clique_solver.cpp)
add_library(pkmc::pkmc ALIAS PKMC)
target_include_directories(PKMC PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_compile_features(PKMC PUBLIC cxx_std_20)

find_package(Threads REQUIRED)
target_link_libraries(PKMC PUBLIC Threads::Threads)
