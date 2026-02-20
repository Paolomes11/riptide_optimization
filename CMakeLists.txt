cmake_minimum_required(VERSION 3.26)
project(riptide_optimization)



set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)



string(APPEND CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Wcast-qual -Wformat=2 -Wundef -Werror=float-equal -Wshadow -Wcast-align -Wunused -Wnull-dereference -Wdouble-promotion -Wimplicit-fallthrough -Wextra-semi -Woverloaded-virtual -Wnon-virtual-dtor -Wold-style-cast")



find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})


add_executable(simulate main.cpp src/physics_list.cpp src/detector_construction.cpp src/action_initialization.cpp src/primary_generator_action.cpp src/sensitive_detector.cpp)
target_link_libraries(simulate PRIVATE ${Geant4_LIBRARIES})
target_include_directories(simulate PRIVATE ${CMAKE_SOURCE_DIR}/include)