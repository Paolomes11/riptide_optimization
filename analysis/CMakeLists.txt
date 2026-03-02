cmake_minimum_required(VERSION 3.26)

# Nome dell'eseguibile per l'analisi
add_executable(plot2D
    plot2D.cpp
)

# Link alle librerie necessarie (ROOT)
find_package(ROOT REQUIRED COMPONENTS
    Core
    RIO
    Hist
    Graf
)
include(${ROOT_USE_FILE})

target_link_libraries(plot2D PRIVATE ${ROOT_LIBRARIES})

# Include directory ROOT se serve
target_include_directories(plot2D PRIVATE ${ROOT_INCLUDE_DIRS})
