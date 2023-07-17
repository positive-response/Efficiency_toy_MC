cmake_minimum_required(VERSION 3.1)

project(Efficiency_toy_MC VERSION 1.0
      LANGUAGES CXX)

add_executable(test.exe test.cpp)
target_compile_options(test.exe PRIVATE -Wall -std=c++14)
