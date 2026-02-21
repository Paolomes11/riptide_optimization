cmake_minimum_required(VERSION 3.26)
project(riptide_optimization)


set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)


string(APPEND CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Wcast-qual -Wformat=2 -Wundef -Werror=float-equal -Wshadow -Wcast-align -Wunused -Wnull-dereference -Wdouble-promotion -Wimplicit-fallthrough -Wextra-semi -Woverloaded-virtual -Wnon-virtual-dtor -Wold-style-cast")


# External dependencies required for the simulation
find_package(XercesC REQUIRED) # XML parsing library for GDML geometry files
find_package(Geant4 REQUIRED ui_all vis_all gdml) # Particle physics simulation toolkit
include(${Geant4_USE_FILE})
find_package(ROOT REQUIRED COMPONENTS Core Tree RIO) # Data analysis framework for ROOT files
## Lyra serve per creare flag personalizzate per eseguire il programma da terminale
# find_package(lyra REQUIRED)           # Command line argument parsing library
# find_package(doctest REQUIRED)        # Lightweight testing framework
find_package(fmt REQUIRED) # Modern C++ formatting library
## spdlog serve per avere log rapidi e personalizzati durante l'esecuzione del programma
find_package(spdlog REQUIRED) # High-performance logging library


add_executable(simulate
    main.cpp
    src/physics_list.cpp
    src/detector_construction.cpp
    src/action_initialization.cpp
    src/primary_generator_action.cpp
    src/sensitive_detector.cpp
    src/run_action.cpp
    src/photocathode_hit.cpp
    src/event_action.cpp
)

target_include_directories(simulate
    PRIVATE
    ${CMAKE_SOURCE_DIR}/include
)

target_link_libraries(simulate
    PRIVATE
    ${Geant4_LIBRARIES}
    fmt::fmt
    spdlog::spdlog
)
