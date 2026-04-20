cmake_minimum_required(VERSION 3.26)

find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Tree)
include(${ROOT_USE_FILE})

add_library(exp_common STATIC
    src/fits_io.cpp
    src/stacking.cpp
)

target_include_directories(exp_common
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE ${ROOT_INCLUDE_DIRS}
)

target_link_libraries(exp_common PUBLIC ${ROOT_LIBRARIES})
target_compile_features(exp_common PUBLIC cxx_std_17)

target_compile_options(exp_common PRIVATE
    -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion
    -Wformat=2 -Wshadow -Wunused -Wnull-dereference -Wold-style-cast
)
