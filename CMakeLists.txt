set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

cmake_minimum_required(VERSION 3.10)
project(hypersonic-reentry-pdtc)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(Matplot++ REQUIRED)
find_package(GTest CONFIG REQUIRED)
enable_testing()

# Clarabel has no vcpkg port (its core is a Rust crate built via Cargo),
# so it is pulled in directly from source and built by its own CMake/Cargo rules.
include(FetchContent)
FetchContent_Declare(
  clarabelcpp
  GIT_REPOSITORY https://github.com/oxfordcontrol/Clarabel.cpp.git
  GIT_TAG        v0.11.1
  GIT_SUBMODULES "Clarabel.rs"
  EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(clarabelcpp)

file(GLOB_RECURSE SOURCE
  "${CMAKE_CURRENT_SOURCE_DIR}/dynamics/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/dynamics/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/spacecraft/*.h"
  "${CMAKE_CURRENT_SOURCE_DIR}/planet/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/planet/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/atmosphere/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/atmosphere/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/attitude/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/attitude/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/aero/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/aero/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/loads/*.h" "${CMAKE_CURRENT_SOURCE_DIR}/loads/*.cpp"
  "${CMAKE_CURRENT_SOURCE_DIR}/main.cpp"
)
# GenerateAeroTable.cpp and GenerateDeflectedGeometry.cpp each have their own
# main() (offline tools, built as separate targets below) -- exclude them
# from the main glob or multiple targets' main()s would collide when linked
# together.
list(REMOVE_ITEM SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/aero/GenerateAeroTable.cpp")
list(REMOVE_ITEM SOURCE "${CMAKE_CURRENT_SOURCE_DIR}/aero/GenerateDeflectedGeometry.cpp")

add_executable(main ${SOURCE})
target_include_directories(main PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/dynamics
  ${CMAKE_CURRENT_SOURCE_DIR}/spacecraft
  ${CMAKE_CURRENT_SOURCE_DIR}/planet
  ${CMAKE_CURRENT_SOURCE_DIR}/atmosphere
  ${CMAKE_CURRENT_SOURCE_DIR}/attitude
  ${CMAKE_CURRENT_SOURCE_DIR}/aero
  ${CMAKE_CURRENT_SOURCE_DIR}/loads
)
target_link_libraries (main Eigen3::Eigen Matplot++::matplot libclarabel_c_shared)

# Offline aero coefficient table generator -- needs only Eigen3, not
# Matplot++/Clarabel. Not run automatically; run it manually to (re)generate
# aero/data/aero_table.csv, which is committed to the repo.
add_executable(aero_table_gen
  aero/GenerateAeroTable.cpp
  aero/PanelGeometry.cpp
  aero/NewtonianAero.cpp
  aero/ShockExpansionAero.cpp
  aero/AeroRegimeDispatch.cpp
  aero/Kriging.cpp
  aero/LatinHypercubeSampler.cpp
  aero/StlMeshLoader.cpp
  aero/SpacecraftGeometry.cpp
)
target_include_directories(aero_table_gen PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/aero
)
target_link_libraries(aero_table_gen PRIVATE Eigen3::Eigen)
# main/unit_tests get C++17 transitively (Matplot++/GTest's CMake targets
# require it), but aero_table_gen only links Eigen3, which doesn't -- needed
# explicitly here for <filesystem> and nested namespace syntax.
target_compile_features(aero_table_gen PRIVATE cxx_std_17)

# Offline deflected-geometry generator -- reads aero/data/doe_points.csv
# (written by aero_table_gen) and rotates each flap's neutral STL about its
# hinge for each DOE point, for handoff to a CFD mesher. Run after
# aero_table_gen has produced doe_points.csv.
add_executable(deflected_geometry_gen
  aero/GenerateDeflectedGeometry.cpp
  aero/PanelGeometry.cpp
  aero/StlMeshLoader.cpp
  aero/NewtonianAero.cpp
  aero/ShockExpansionAero.cpp
  aero/AeroRegimeDispatch.cpp
)
target_include_directories(deflected_geometry_gen PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/aero
)
target_link_libraries(deflected_geometry_gen PRIVATE Eigen3::Eigen)
target_compile_features(deflected_geometry_gen PRIVATE cxx_std_17)

add_executable(unit_tests
  tests/AtmosphereModelTests.cpp
  tests/DynamicsTests.cpp
  tests/AeroTests.cpp
  tests/HeatingLoadModelTests.cpp
  atmosphere/AtmosphereModel.cpp
  dynamics/DescentDynamics.cpp
  attitude/QuaternionUtils.cpp
  planet/PlanetConfig.cpp
  aero/PanelGeometry.cpp
  aero/NewtonianAero.cpp
  aero/ShockExpansionAero.cpp
  aero/AeroRegimeDispatch.cpp
  aero/Kriging.cpp
  aero/LatinHypercubeSampler.cpp
  aero/StlMeshLoader.cpp
  aero/SpacecraftGeometry.cpp
  aero/AeroAngles.cpp
  aero/AeroCoefficientTable.cpp
  loads/HeatingLoadModel.cpp
)
target_include_directories(unit_tests PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}/dynamics
  ${CMAKE_CURRENT_SOURCE_DIR}/spacecraft
  ${CMAKE_CURRENT_SOURCE_DIR}/planet
  ${CMAKE_CURRENT_SOURCE_DIR}/atmosphere
  ${CMAKE_CURRENT_SOURCE_DIR}/attitude
  ${CMAKE_CURRENT_SOURCE_DIR}/aero
  ${CMAKE_CURRENT_SOURCE_DIR}/loads
)
target_link_libraries(unit_tests PRIVATE Eigen3::Eigen GTest::gtest GTest::gtest_main)
target_compile_definitions(unit_tests PRIVATE TESTS_SOURCE_DIR="${CMAKE_CURRENT_SOURCE_DIR}")
include(GoogleTest)
gtest_discover_tests(unit_tests)

if(WIN32)
  # Clarabel.rs's CMake only distinguishes Debug/Release via CMAKE_BUILD_TYPE,
  # which multi-config generators (e.g. Visual Studio) never set, so it always
  # builds/outputs to rust_wrapper/target/debug here.
  add_custom_command(
    TARGET main POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy_if_different
      "${clarabelcpp_SOURCE_DIR}/rust_wrapper/target/debug/clarabel_c.dll"
      "$<TARGET_FILE_DIR:main>"
  )
endif()