cmake_minimum_required(VERSION 3.26)

# ROOT
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Tree)
include(${ROOT_USE_FILE})

# Libreria statica exp1_lib
add_library(exp1_lib STATIC
    src/fit_reader.cpp
    src/frame_stats.cpp
    src/image_analysis.cpp
)

target_include_directories(exp1_lib
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE ${ROOT_INCLUDE_DIRS}
)

target_link_libraries(exp1_lib
    PUBLIC ${ROOT_LIBRARIES}
)

target_compile_features(exp1_lib PUBLIC cxx_std_17)

set(COMMON_WARNINGS
    -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion
    -Wformat=2 -Wshadow -Wunused -Wnull-dereference
    -Wold-style-cast
)
target_compile_options(exp1_lib PRIVATE ${COMMON_WARNINGS})

# Eseguibile principale
add_executable(exp1_main exp1_main.cpp)

# Eseguibile per visualizzare il .root
add_executable(exp1_root_viewer exp1_root_viewer.cpp)

target_link_libraries(exp1_main
    PRIVATE exp1_lib ${ROOT_LIBRARIES}
)

target_link_libraries(exp1_root_viewer
    PRIVATE ${ROOT_LIBRARIES}
)

target_include_directories(exp1_main
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${ROOT_INCLUDE_DIRS}
)

target_include_directories(exp1_root_viewer
    PRIVATE ${ROOT_INCLUDE_DIRS}
)

target_compile_features(exp1_main PRIVATE cxx_std_17)
target_compile_features(exp1_root_viewer PRIVATE cxx_std_17)
target_compile_options(exp1_main PRIVATE ${COMMON_WARNINGS})
target_compile_options(exp1_root_viewer PRIVATE ${COMMON_WARNINGS})

# Directory di output
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/output/exp1)

# Messaggio informativo
message(STATUS "exp1 analysis module configured")
message(STATUS "  Place FITS files in: analysis/exp1/data/{good,bad,background}/")
message(STATUS "  Or on SSD:           <ssd-mount>/exp1/{good,bad,background}/")
