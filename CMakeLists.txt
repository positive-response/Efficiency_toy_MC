cmake_minimum_required(VERSION 3.1)

project(Efficiency_toy_MC VERSION 1.0
      LANGUAGES CXX)

add_executable(test.exe test.cpp)
target_compile_options(test.exe PRIVATE -Wall -std=c++14)

find_package(ROOT 6.18 CONFIG REQUIRED COMPONENTS Core RIO Hist RooStats Graf Graf3d ROOTDataFrame)
include("${ROOT_USE_FILE}")


add_executable(phaspace4.exe phasespace4.cpp)
target_link_libraries(phaspace4.exe PUBLIC ROOT::Physics ROOT::Core ROOT::Tree ROOT::RIO ROOT::Net ROOT::Hist ROOT::RooStats ROOT::Graf ROOT::Graf3d ROOT::ROOTDataFrame)
target_compile_options(phaspace4.exe PRIVATE -Wall -std=c++14)

add_executable(klein.exe Klein_4.cpp)
target_link_libraries(klein.exe PUBLIC ROOT::Physics ROOT::Core ROOT::Tree ROOT::RIO ROOT::Net ROOT::Hist ROOT::RooStats ROOT::Graf ROOT::Graf3d ROOT::ROOTDataFrame)
target_compile_options(klein.exe PRIVATE -Wall -std=c++14)

add_executable(inter.exe inter.cpp)
target_compile_options(inter.exe PRIVATE -Wall -std=c++14)
target_link_libraries(inter.exe PUBLIC ROOT::Physics ROOT::Core ROOT::Tree ROOT::RIO ROOT::Net ROOT::Hist ROOT::RooStats ROOT::Graf ROOT::Graf3d ROOT::ROOTDataFrame)


