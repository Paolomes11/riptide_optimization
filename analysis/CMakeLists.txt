cmake_minimum_required(VERSION 3.26)

# Eseguibile esistente per l'analisi 2D
add_executable(plot2D
    plot2D.cpp
)

# Nuovo eseguibile per l'analisi 3D
add_executable(plot3D
    plot3D.cpp
)

# Eseguibile per il beam scan
add_executable(beam_scan_plot
    beam_scan_plot.cpp
)

# Eseguibile per la mappa di covarianza e media
add_executable(m_c_creator
    m_c_creator.cpp
)

# Eseguibile per l'estrazione di media e matrice di covarianza PSF per tutti i run
add_executable(psf_extractor
    psf_extractor.cpp
)

# ROOT
find_package(ROOT REQUIRED COMPONENTS
    Core
    RIO
    Hist
    Graf
    Graf3d
)
include(${ROOT_USE_FILE})

# Link delle librerie ROOT agli eseguibili
target_link_libraries(plot2D PRIVATE ${ROOT_LIBRARIES})
target_link_libraries(plot3D PRIVATE ${ROOT_LIBRARIES})
target_link_libraries(beam_scan_plot PRIVATE ${ROOT_LIBRARIES})
target_link_libraries(m_c_creator PRIVATE ${ROOT_LIBRARIES})
target_link_libraries(psf_extractor PRIVATE ${ROOT_LIBRARIES})

# Include directory ROOT se serve
target_include_directories(plot2D PRIVATE ${ROOT_INCLUDE_DIRS})
target_include_directories(plot3D PRIVATE ${ROOT_INCLUDE_DIRS})
target_include_directories(beam_scan_plot PRIVATE ${ROOT_INCLUDE_DIRS})
target_include_directories(m_c_creator PRIVATE ${ROOT_INCLUDE_DIRS})
target_include_directories(psf_extractor PRIVATE ${ROOT_INCLUDE_DIRS})

# Sottocartella per analisi PSF
add_subdirectory(psf_analysis)

# Sottocartella per esperimento 1 (immagini FITS)
add_subdirectory(exp1)
