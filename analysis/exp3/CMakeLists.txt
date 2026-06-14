cmake_minimum_required(VERSION 3.26)

find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf Tree)
include(${ROOT_USE_FILE})

find_path(JSON_INCLUDE_DIR nlohmann/json.hpp
    HINTS ${ROOT_INCLUDE_DIRS} /usr/include /usr/local/include REQUIRED)

add_library(exp3_lib STATIC
    src/exp3_analysis.cpp
    src/exp3b_analysis.cpp
    src/exp3_plots.cpp
    src/exp3_config.cpp
)

target_include_directories(exp3_lib
    PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}/include
    PRIVATE ${ROOT_INCLUDE_DIRS}
            ${JSON_INCLUDE_DIR}
            ${CMAKE_SOURCE_DIR}/analysis/psf_analysis
            ${CMAKE_SOURCE_DIR}/analysis/exp2/include
            ${CMAKE_SOURCE_DIR}/analysis/exp_common/include
)

# exp2_lib porta già PUBLIC: exp_common, psf_analysis, ROOT_LIBRARIES
target_link_libraries(exp3_lib
    PUBLIC exp2_lib
           ${ROOT_LIBRARIES}
)

target_compile_features(exp3_lib PUBLIC cxx_std_17)

set(COMMON_WARNINGS
    -Wall -Wextra -Wpedantic -Wconversion -Wsign-conversion
    -Wformat=2 -Wshadow -Wunused -Wnull-dereference -Wold-style-cast
)
target_compile_options(exp3_lib PRIVATE ${COMMON_WARNINGS})

add_executable(exp3_main exp3_main.cpp)
target_link_libraries(exp3_main PRIVATE exp3_lib ${ROOT_LIBRARIES})
target_include_directories(exp3_main
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include
            ${ROOT_INCLUDE_DIRS}
)
target_compile_features(exp3_main PRIVATE cxx_std_17)
target_compile_options(exp3_main PRIVATE ${COMMON_WARNINGS})

file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/output/exp3)
file(MAKE_DIRECTORY ${CMAKE_SOURCE_DIR}/output/exp3/calib)

# Unit test: exp3 calibrazione + tracce parallele + exp3b
add_executable(test_exp3_calib test_exp3_calib.cpp)
target_link_libraries(test_exp3_calib
    PRIVATE exp3_lib
            ${ROOT_LIBRARIES})
target_include_directories(test_exp3_calib
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include
            ${ROOT_INCLUDE_DIRS}
            ${JSON_INCLUDE_DIR})
target_compile_features(test_exp3_calib PRIVATE cxx_std_17)
target_compile_options(test_exp3_calib  PRIVATE ${COMMON_WARNINGS})

enable_testing()
add_test(NAME exp3_calib_test
         COMMAND test_exp3_calib
         WORKING_DIRECTORY ${CMAKE_SOURCE_DIR})
set_tests_properties(exp3_calib_test PROPERTIES
    LABELS "unit;exp3"
    TIMEOUT 30)

message(STATUS "exp3 analysis module configured (v2: parallel traces + line calib)")
message(STATUS "  Place FITS files in: data/exp3/{calib,good,bad}/d{dist}/")
