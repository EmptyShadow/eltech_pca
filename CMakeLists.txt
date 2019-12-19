cmake_minimum_required(VERSION 3.15)
project(pca)
set(CMAKE_CXX_STANDARD 14)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

find_package(MPI REQUIRED)

include_directories(SYSTEM ${MPI_INCLUDE_PATH})

add_executable(lab1
        src/lab1/main.cpp)
target_link_libraries(lab1 PUBLIC OpenMP::OpenMP_CXX)

add_executable(lab2
        src/lab2/main.cpp)
target_link_libraries(lab2 PUBLIC OpenMP::OpenMP_CXX)

add_executable(lab2_mpich
        src/lab2_mpich/main.cpp)
target_link_libraries(lab2_mpich ${MPI_C_LIBRARIES})

if(MPI_COMPILE_FLAGS)
    set_target_properties(lab2_mpich PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(lab2_mpich PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

add_executable(lab3
        src/lab3/main.cpp)
target_link_libraries(lab3 PUBLIC OpenMP::OpenMP_CXX)