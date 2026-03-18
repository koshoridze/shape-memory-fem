cmake_minimum_required(VERSION 3.16)

project(shape_memory_fem LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED)

add_executable(shape_memory_fem
    heat_solver.cpp
    importmesh.cpp
    linelastity_solver.cpp
    main.cpp
    shapefun.cpp
    sma_solver.cpp
    solver.cpp
    vtkexport.cpp
)

target_include_directories(shape_memory_fem PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}
)

target_link_libraries(shape_memory_fem PRIVATE
    Eigen3::Eigen
)
