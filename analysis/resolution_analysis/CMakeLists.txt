cmake_minimum_required(VERSION 3.26)

add_executable(resolution_map
    resolution_map.cpp
)

target_link_libraries(resolution_map PRIVATE ${ROOT_LIBRARIES})
target_include_directories(resolution_map PRIVATE ${ROOT_INCLUDE_DIRS})

