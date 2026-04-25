cmake_minimum_required(VERSION 3.26)

find_package(ROOT REQUIRED COMPONENTS Core RIO Hist Graf)
include(${ROOT_USE_FILE})

# Eseguibile principale
add_executable(pareto_selector
    pareto_selector.cpp
)
target_link_libraries(pareto_selector
    PRIVATE ${ROOT_LIBRARIES} spdlog::spdlog
)
target_include_directories(pareto_selector
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
            ${ROOT_INCLUDE_DIRS}
)
target_compile_features(pareto_selector PRIVATE cxx_std_17)

# Test unitari (pareto_core.hpp pura, nessun ROOT richiesto a runtime)
add_executable(test_pareto_selector
    test_pareto_selector.cpp
)
target_link_libraries(test_pareto_selector
    PRIVATE ${ROOT_LIBRARIES}
)
target_include_directories(test_pareto_selector
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}
            ${ROOT_INCLUDE_DIRS}
)
target_compile_features(test_pareto_selector PRIVATE cxx_std_17)

enable_testing()
add_test(
    NAME pareto_selector_unit
    COMMAND test_pareto_selector
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
)
set_tests_properties(pareto_selector_unit PROPERTIES
    LABELS "unit;pareto_analysis"
    TIMEOUT 30
)
