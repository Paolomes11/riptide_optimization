cmake_minimum_required(VERSION 3.26)

find_package(OpenMP REQUIRED)

add_executable(dof_plot
    dof_plot.cpp
)

add_executable(dof_map
    dof_map.cpp
)

add_executable(magnification_map
    magnification_map.cpp
)

add_executable(psf_dof_extractor
    psf_dof_extractor.cpp
)

target_link_libraries(dof_plot PRIVATE ${ROOT_LIBRARIES})
target_link_libraries(dof_map PRIVATE ${ROOT_LIBRARIES} OpenMP::OpenMP_CXX)
target_link_libraries(magnification_map PRIVATE ${ROOT_LIBRARIES})
target_link_libraries(psf_dof_extractor PRIVATE ${ROOT_LIBRARIES} OpenMP::OpenMP_CXX)

target_include_directories(dof_plot PRIVATE ${ROOT_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_include_directories(dof_map PRIVATE ${ROOT_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_include_directories(magnification_map PRIVATE ${ROOT_INCLUDE_DIRS} ${CMAKE_CURRENT_SOURCE_DIR}/..)
target_include_directories(psf_dof_extractor PRIVATE ${ROOT_INCLUDE_DIRS})
