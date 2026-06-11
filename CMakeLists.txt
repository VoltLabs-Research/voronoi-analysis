cmake_minimum_required(VERSION 3.20 FATAL_ERROR)

project(voronoi-analysis VERSION 1.0.2 LANGUAGES CXX)

find_package(coretoolkit REQUIRED)

volt_add_plugin(voronoi-analysis
    MAIN_SOURCE src/plugin.cpp
)
