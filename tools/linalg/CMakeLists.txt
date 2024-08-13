cmake_minimum_required(VERSION 3.10)

project(LinAlg)
set(CMAKE_CXX_STANDARD 20)
set(PYBIND11_PYTHON_VERSION 3.6)
set(CMAKE_CXX_FLAGS "-Wall -Wextra -fPIC")
list(APPEND CMAKE_PREFIX_PATH "/home/40265864@ecit.qub.ac.uk/.local/lib/python3.10/site-packages/pybind11")

find_package(pybind11 REQUIRED)
find_package(Eigen3 REQUIRED)
find_package(OpenMP)

pybind11_add_module(${PROJECT_NAME} ../linalgtools.cpp)

target_compile_definitions(${PROJECT_NAME} PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
target_include_directories(${PROJECT_NAME} PRIVATE ${PYBIND11_INCLUDE_DIRS})
target_link_libraries(${PROJECT_NAME} PRIVATE Eigen3::Eigen)
target_link_libraries(${PROJECT_NAME}  PUBLIC OpenMP::OpenMP_CXX)
