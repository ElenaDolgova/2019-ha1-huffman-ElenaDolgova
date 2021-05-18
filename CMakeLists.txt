cmake_minimum_required(VERSION 3.13)
project(haffman)

set(CMAKE_CXX_STANDARD 14)

add_executable(huffman main.cpp huffman.cpp huffman.hpp)