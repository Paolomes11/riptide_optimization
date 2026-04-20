cmake_minimum_required(VERSION 3.26)

find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Tree)
include(${ROOT_USE_FILE})

add_library(exp2_lib STATIC
    src/trace_extractor.cpp
    src/exp2_analysis.cpp
)

target_include_directories(exp2_lib
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE ${ROOT_INCLUDE_DIRS}
    ${CMAKE_SOURCE_DIR}/analysis/psf_analysis
)

target_link_libraries(exp2_lib
    PUBLIC exp_common
    psf_analysis
    ${ROOT_LIBRARIES}
)

target_compile_features(exp2_lib PUBLIC cxx_std_17)

set(COMMON_WARNINGS
    -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion
    -Wformat=2 -Wshadow -Wunused -Wnull-dereference -Wold-style-cast
)
target_compile_options(exp2_lib PRIVATE ${COMMON_WARNINGS})

add_executable(exp2_main exp2_main.cpp)
target_link_libraries(exp2_main PRIVATE exp2_lib ${ROOT_LIBRARIES})
target_include_directories(exp2_main PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${ROOT_INCLUDE_DIRS}
)
target_compile_features(exp2_main PRIVATE cxx_std_17)
target_compile_options(exp2_main PRIVATE ${COMMON_WARNINGS})

file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/output/exp2)

message(STATUS "exp2 analysis module configured")
message(STATUS "  Data: analysis/exp2/data/{good_focus,good_nofocus,bad_focus,bad_nofocus,background}/")
