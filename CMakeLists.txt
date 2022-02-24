cmake_minimum_required(VERSION 3.10)
project(ml-predictions)

find_library(TENSORFLOW_LIB tensorflow HINT ./lib/libtensorflow2/lib)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer -fsanitize=address")
set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -lasan")

add_executable(ml-predictions ml-predictions.cc)
target_include_directories(ml-predictions PRIVATE ../../include ./lib/libtensorflow2/include)
target_link_libraries (ml-predictions "${TENSORFLOW_LIB}")
