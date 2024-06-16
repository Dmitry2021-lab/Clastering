cmake_minimum_required(VERSION 3.24)
project(Clustering)

set(CMAKE_CXX_STANDARD 14)

add_executable(Clustering MainPy.cpp GaussianMixture.cpp GaussianMixture.h LindeBuzoGray.cpp LindeBuzoGray.h pybind11.h)
