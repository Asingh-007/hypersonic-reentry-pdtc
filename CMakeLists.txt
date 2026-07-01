set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

cmake_minimum_required(VERSION 3.10)
project(hypersonic-reentry-pdtc)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

find_package(Eigen3 REQUIRED NO_MODULE)

file(GLOB SOURCE 
"*.h"
"*.cpp"
)
add_executable(main ${SOURCE})
target_link_libraries (main Eigen3::Eigen)