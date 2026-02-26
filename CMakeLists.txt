cmake_minimum_required(VERSION 3.26)
project(riptide_optimization)



set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)



# Geant4
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

# ROOT
find_package(ROOT REQUIRED COMPONENTS
    Core
    Hist
    RIO
)

include(${ROOT_USE_FILE})

# spdlog
find_package(spdlog REQUIRED)

# lyra
add_library(lyra INTERFACE)
target_include_directories(lyra INTERFACE
    ${CMAKE_SOURCE_DIR}/external
)

# Eseguibile principale
add_executable(simulate
    main.cpp
    src/physics_list.cpp
    src/detector_construction.cpp
    src/action_initialization.cpp
    src/primary_generator_action.cpp
    src/sensitive_detector.cpp
    src/run_action.cpp
    src/event_action.cpp
    src/optimizer.cpp
    src/efficiency_collector.cpp
)

# Warnings specifici per simulate
target_compile_options(simulate PRIVATE
    -Wall
    -Wextra
    -Wpedantic
    -Wconversion
    -Wsign-conversion
    -Wcast-qual
    -Wformat=2
    -Wundef
    -Werror=float-equal
    -Wshadow
    -Wcast-align
    -Wunused
    -Wnull-dereference
    -Wdouble-promotion
    -Wimplicit-fallthrough
    -Wextra-semi
    -Woverloaded-virtual
    -Wnon-virtual-dtor
    -Wold-style-cast
)

# Include directories
target_include_directories(simulate PRIVATE ${CMAKE_SOURCE_DIR}/include)

# Link alle librerie
target_link_libraries(simulate PRIVATE
    ${Geant4_LIBRARIES}
    ${ROOT_LIBRARIES}
    spdlog::spdlog
    lyra
)

# Crea la cartella di output
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/output)

# Subdirectory per l'analisi dei dati (debug e sviluppo)
# add_subdirectory(analysis)
