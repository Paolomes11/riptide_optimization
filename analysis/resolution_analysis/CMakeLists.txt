cmake_minimum_required(VERSION 3.26)

find_package(OpenMP REQUIRED)

add_executable(resolution_map
    resolution_map.cpp
)

target_link_libraries(resolution_map PRIVATE ${ROOT_LIBRARIES} OpenMP::OpenMP_CXX)
target_include_directories(resolution_map PRIVATE ${ROOT_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/..)

