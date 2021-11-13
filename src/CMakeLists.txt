cmake_minimum_required(VERSION 3.10)

project(Laplace)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

add_executable(Laplace_DirichletBC TestCase1.cpp SolveGaussSeidel.h SimulationParameter.h VtkData.h VtkData.cpp)

add_executable(Laplace_NeumannBC TestCase2.cpp SolveGaussSeidel.h SimulationParameter.h VtkData.h VtkData.cpp)