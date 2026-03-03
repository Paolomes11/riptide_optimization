cmake_minimum_required(VERSION 3.26)
project(riptide_optimization)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# -----------------------------
# Librerie esterne
# -----------------------------
# Geant4
find_package(Geant4 REQUIRED ui_all vis_all)
include(${Geant4_USE_FILE})

# ROOT
find_package(ROOT REQUIRED COMPONENTS Core Hist RIO)
include(${ROOT_USE_FILE})

# spdlog
find_package(spdlog REQUIRED)

# lyra
add_library(lyra INTERFACE)
target_include_directories(lyra INTERFACE ${CMAKE_SOURCE_DIR}/external)

# -----------------------------
# Libreria comune
# -----------------------------
file(GLOB COMMON_SRC src/common/*.cpp)
add_library(common_lib ${COMMON_SRC})
target_include_directories(common_lib PUBLIC ${CMAKE_SOURCE_DIR}/include/common)
target_link_libraries(common_lib PUBLIC ${ROOT_LIBRARIES} spdlog::spdlog lyra ${Geant4_LIBRARIES})

# -----------------------------
# Directory include globali
# -----------------------------
include_directories(
    ${CMAKE_SOURCE_DIR}/include
    ${CMAKE_SOURCE_DIR}/include/common
)

# -----------------------------
# Warnings comuni
# -----------------------------
set(COMMON_WARNINGS
    -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion -Wcast-qual
    -Wformat=2 -Wundef -Werror=float-equal -Wshadow -Wcast-align
    -Wunused -Wnull-dereference -Wdouble-promotion -Wimplicit-fallthrough
    -Wextra-semi -Woverloaded-virtual -Wnon-virtual-dtor -Wold-style-cast
)

# -----------------------------
# Costruzione dei progetti (main)
# -----------------------------
file(GLOB PROGRAM_MAIN programs/*_main.cpp)

foreach(MAIN_FILE ${PROGRAM_MAIN})
    get_filename_component(PROG_NAME ${MAIN_FILE} NAME_WE)

    # Nome sottocartella sorgenti della simulazione
    string(REPLACE "_main" "" PROG_SUBDIR ${PROG_NAME})
    set(SRC_DIR src/${PROG_SUBDIR})

    # Libreria specifica della simulazione
    file(GLOB SPEC_SRC ${SRC_DIR}/*.cpp)
    if(SPEC_SRC)
        add_library(${PROG_NAME}_lib ${SPEC_SRC})
        target_include_directories(${PROG_NAME}_lib PUBLIC ${CMAKE_SOURCE_DIR}/include/${PROG_SUBDIR})
        target_link_libraries(${PROG_NAME}_lib PUBLIC common_lib)
        target_compile_options(${PROG_NAME}_lib PRIVATE ${COMMON_WARNINGS})
    endif()

    # Eseguibile
    add_executable(${PROG_NAME} ${MAIN_FILE})
    target_include_directories(${PROG_NAME} PRIVATE ${CMAKE_SOURCE_DIR}/include)
    target_link_libraries(${PROG_NAME} PRIVATE common_lib ${PROG_NAME}_lib)
    target_compile_options(${PROG_NAME} PRIVATE ${COMMON_WARNINGS})
endforeach()

# -----------------------------
# Crea la cartella di output
# -----------------------------
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/output)

# -----------------------------
# Subdirectory per l'analisi dei dati
# -----------------------------
add_subdirectory(analysis)