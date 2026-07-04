set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

cmake_minimum_required(VERSION 3.10)
project(hypersonic-reentry-pdtc)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

find_package(Eigen3 REQUIRED NO_MODULE)
find_package(Matplot++ REQUIRED)

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

file(GLOB SOURCE
"*.h"
"*.cpp"
)
add_executable(main ${SOURCE})
target_link_libraries (main Eigen3::Eigen Matplot++::matplot libclarabel_c_shared)

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