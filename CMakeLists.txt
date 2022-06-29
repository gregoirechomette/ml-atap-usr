cmake_minimum_required(VERSION 3.10)
project(main)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fno-omit-frame-pointer")
set(CMAKE_LINKER_FLAGS "${CMAKE_LINKER_FLAGS} -lasan")

add_library(affected-population STATIC model.cc popGrid.cc)

set(APP_EXE main-test)

add_executable(${APP_EXE} main.cc)

target_link_libraries(${APP_EXE} affected-population)
