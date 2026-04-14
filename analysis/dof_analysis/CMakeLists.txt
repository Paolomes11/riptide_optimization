cmake_minimum_required(VERSION 3.26)

add_executable(dof_plot
    dof_plot.cpp
)

add_executable(dof_map
    dof_map.cpp
)

target_link_libraries(dof_plot PRIVATE ${ROOT_LIBRARIES})
target_link_libraries(dof_map PRIVATE ${ROOT_LIBRARIES})

target_include_directories(dof_plot PRIVATE ${ROOT_INCLUDE_DIRS})
target_include_directories(dof_map PRIVATE ${ROOT_INCLUDE_DIRS})
