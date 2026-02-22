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
