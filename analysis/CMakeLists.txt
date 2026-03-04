cmake_minimum_required(VERSION 3.26)

# Eseguibile esistente per l'analisi 2D
add_executable(plot2D
    plot2D.cpp
)

# Nuovo eseguibile per il beam scan
add_executable(beam_scan_plot
    beam_scan_plot.cpp
)

# Nuovo eseguibile per il beam scan
add_executable(geant4_heatmap
    geant4_heatmap.cpp
)


# ROOT
find_package(ROOT REQUIRED COMPONENTS
    Core
    RIO
    Hist
    Graf
)
include(${ROOT_USE_FILE})

# Link delle librerie ROOT agli eseguibili
target_link_libraries(plot2D PRIVATE ${ROOT_LIBRARIES})
target_link_libraries(beam_scan_plot PRIVATE ${ROOT_LIBRARIES})
target_link_libraries(geant4_heatmap PRIVATE ${ROOT_LIBRARIES})

# Include directory ROOT se serve
target_include_directories(plot2D PRIVATE ${ROOT_INCLUDE_DIRS})
target_include_directories(beam_scan_plot PRIVATE ${ROOT_INCLUDE_DIRS})
target_include_directories(geant4_heatmap PRIVATE ${ROOT_INCLUDE_DIRS})