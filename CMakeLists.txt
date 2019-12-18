cmake_minimum_required(VERSION 3.15)
project(pca)
set(CMAKE_CXX_STANDARD 14)

find_package(MPI REQUIRED)
find_package(OpenMP)

include_directories(${MPI_INCLUDE_PATH})

add_executable(lab1
        src/lab1/main.cpp)
target_link_libraries(lab1 PUBLIC OpenMP::OpenMP_CXX)

add_executable(lab2
        src/lab2/main.cpp)
target_link_libraries(lab2 PUBLIC OpenMP::OpenMP_CXX ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
    set_target_properties(lab2 PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(lab2 PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()