cmake_minimum_required(VERSION 3.10)
project(ml-predictions)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer")
set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -lasan")

add_executable(ml-predictions ml-predictions.cc input.cc model.cc popGrid.cc)