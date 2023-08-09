cmake_minimum_required(VERSION 3.1)

project(Efficiency_toy_MC VERSION 1.0
      LANGUAGES CXX)

add_executable(test.exe test.cpp)
target_compile_options(test.exe PRIVATE -Wall -std=c++14)

find_package(ROOT 6.18 CONFIG REQUIRED COMPONENTS Core RIO Hist RooStats Graf Graf3d ROOTDataFrame)
include("${ROOT_USE_FILE}")

add_executable(GeometricalEfficiency.exe GeometricalEfficiency.cpp)
target_link_libraries(GeometricalEfficiency.exe PUBLIC ROOT::Physics ROOT::Core ROOT::Tree ROOT::RIO ROOT::Net ROOT::Hist ROOT::RooStats ROOT::Graf ROOT::Graf3d ROOT::ROOTDataFrame)
target_compile_options(GeometricalEfficiency.exe PRIVATE -Wall -std=c++14)

add_executable(RegistrationEfficiency.exe RegistrationEfficiency.cpp )
target_link_libraries(RegistrationEfficiency.exe PUBLIC ROOT::Physics ROOT::Core ROOT::Tree ROOT::RIO ROOT::Net ROOT::Hist ROOT::RooStats ROOT::Graf ROOT::Graf3d ROOT::ROOTDataFrame)
target_compile_options(RegistrationEfficiency.exe PRIVATE -Wall -std=c++14)

##add_executable(DetectionEfficiency.exe DetectionEfficiency.cpp RegistrationEfficiency.h)
##target_compile_options(DetectionEfficiency.exe PRIVATE -Wall -std=c++14)
##target_link_libraries(DetectionEfficiency.exe PUBLIC ROOT::Physics ROOT::Core ROOT::Tree ROOT::RIO ROOT::Net ROOT::Hist ROOT::RooStats ROOT::Graf ROOT::Graf3d ROOT::ROOTDataFrame)


