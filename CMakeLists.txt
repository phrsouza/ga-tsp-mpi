cmake_minimum_required(VERSION 3.6)
project(TP1_nqueens)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})


#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -pthread")
set(CMAKE_C_COMPILER mpicc)

set(SOURCE_FILES main.c)
add_executable(PPD ${SOURCE_FILES})