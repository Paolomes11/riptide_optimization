cmake_minimum_required(VERSION 3.16)
project(riptide_analysis)

find_package(ROOT REQUIRED COMPONENTS
    Core
    RIO
    Tree
    Hist
    Graf
    Gpad
)

add_executable(viewer
    src/viewer.cpp
)

target_link_libraries(viewer
    ROOT::Core
    ROOT::RIO
    ROOT::Tree
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad
)

# Eseguibile per efficienza geometrica
add_executable(geometric_efficiency
    src/efficiency.cpp
    main_efficiency.cpp
)

target_include_directories(geometric_efficiency PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)

target_link_libraries(geometric_efficiency
    ROOT::Core
    ROOT::RIO
    ROOT::Tree
    ROOT::Hist
    ROOT::Graf
    ROOT::Gpad
)
