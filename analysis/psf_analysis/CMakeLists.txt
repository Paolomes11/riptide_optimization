cmake_minimum_required(VERSION 3.26)

# ROOT
find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf)
include(${ROOT_USE_FILE})

# nlohmann/json (header-only; bundled con ROOT o installato separatamente)
find_path(JSON_INCLUDE_DIR nlohmann/json.hpp
    HINTS ${ROOT_INCLUDE_DIRS} /usr/include /usr/local/include
    REQUIRED
)

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

# Eseguibile q_map
add_executable(q_map
    q_map.cpp
)
target_link_libraries(q_map
    PRIVATE psf_analysis ${ROOT_LIBRARIES} spdlog::spdlog
)
target_include_directories(q_map
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    ${ROOT_INCLUDE_DIRS}
    ${JSON_INCLUDE_DIR}
)
target_compile_features(q_map PRIVATE cxx_std_17)

# Eseguibile chi2_map
add_executable(chi2_map
    chi2_map.cpp
)
target_link_libraries(chi2_map
    PRIVATE psf_analysis ${ROOT_LIBRARIES} spdlog::spdlog
)
target_include_directories(chi2_map
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    ${ROOT_INCLUDE_DIRS}
    ${JSON_INCLUDE_DIR}
)
target_compile_features(chi2_map PRIVATE cxx_std_17)

# Test unitari fit_trace e compute_Q
add_executable(test_fit_trace
    test_fit_trace.cpp
)
target_link_libraries(test_fit_trace
    PRIVATE psf_analysis
    ${ROOT_LIBRARIES}
)
target_include_directories(test_fit_trace
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    ${ROOT_INCLUDE_DIRS}
)
target_compile_features(test_fit_trace PRIVATE cxx_std_17)

# Registra il test in CTest in modo che `ctest` lo esegua automaticamente
# dalla directory di build del progetto principale.
enable_testing()
add_test(
    NAME fit_trace_unit
    COMMAND test_fit_trace
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
# Il test e' considerato passato se il processo esce con codice 0.
set_tests_properties(fit_trace_unit PROPERTIES
    LABELS "unit;psf_analysis"
    TIMEOUT 30
)

add_executable(test_chi2_diagnostics
    test_chi2_diagnostics.cpp
)
target_link_libraries(test_chi2_diagnostics
    PRIVATE psf_analysis
    ${ROOT_LIBRARIES}
)
target_include_directories(test_chi2_diagnostics
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
    ${ROOT_INCLUDE_DIRS}
)
target_compile_features(test_chi2_diagnostics PRIVATE cxx_std_17)

add_test(
    NAME chi2_diagnostics
    COMMAND test_chi2_diagnostics
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
set_tests_properties(chi2_diagnostics PROPERTIES
    LABELS "unit;psf_analysis;chi2"
    TIMEOUT 30
)
