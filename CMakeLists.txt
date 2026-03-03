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

# Libreria comune
file(GLOB_RECURSE COMMON_SRC src/common/*.cpp)
add_library(common_lib ${COMMON_SRC})
target_include_directories(common_lib PUBLIC ${CMAKE_SOURCE_DIR}/include/common)
target_link_libraries(common_lib PUBLIC ${ROOT_LIBRARIES} spdlog::spdlog lyra ${Geant4_LIBRARIES})

# Rileva automaticamente tutti i main in programs/
file(GLOB PROGRAM_MAIN programs/*_main.cpp)

# Directory include
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include/common
    ${CMAKE_SOURCE_DIR}/include/optimization
)

foreach(MAIN_FILE ${PROGRAM_MAIN})
    get_filename_component(PROG_NAME ${MAIN_FILE} NAME_WE)
    # Trova eventuali src specifici della simulazione
    string(REPLACE "_main" "" PROG_SUBDIR ${PROG_NAME})
    file(GLOB_RECURSE SPEC_SRC src/${PROG_SUBDIR}/*.cpp)

    add_executable(${PROG_NAME} ${MAIN_FILE} ${SPEC_SRC})
    target_include_directories(${PROG_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/include)
    target_link_libraries(${PROG_NAME} PRIVATE common_lib)
    
    # Warnings specifici (puoi aggiungerli anche qui)
    target_compile_options(${PROG_NAME} PRIVATE
        -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Wcast-qual
        -Wformat=2 -Wundef -Werror=float-equal -Wshadow -Wcast-align
        -Wunused -Wnull-dereference -Wdouble-promotion -Wimplicit-fallthrough
        -Wextra-semi -Woverloaded-virtual -Wnon-virtual-dtor -Wold-style-cast
    )
endforeach()


# Crea la cartella di output
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/output)

# Subdirectory per l'analisi dei dati (debug e sviluppo)
add_subdirectory(analysis)
