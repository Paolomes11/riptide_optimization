cmake_minimum_required(VERSION 3.26)

# ROOT
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf)
include(${ROOT_USE_FILE})

# Libreria psf_analysis
add_library(psf_analysis STATIC
    psf_interpolator.cpp
)
target_include_directories(psf_analysis
    PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
    PRIVATE ${ROOT_INCLUDE_DIRS}
)
target_link_libraries(psf_analysis PUBLIC ${ROOT_LIBRARIES})
target_compile_features(psf_analysis PUBLIC cxx_std_17)

# Eseguibile trace_viewer
add_executable(trace_viewer
    trace_viewer.cpp
)
target_link_libraries(trace_viewer
    PRIVATE psf_analysis
    ${ROOT_LIBRARIES}
)
target_include_directories(trace_viewer
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    ${ROOT_INCLUDE_DIRS}
)
target_compile_features(trace_viewer PRIVATE cxx_std_17)
