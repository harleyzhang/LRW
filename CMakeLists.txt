#This project implement the read and write gdf file
project(clustering)
cmake_minimum_required(VERSION 2.8)

#make sure that gcc is used by cmake. This can be done by 
#setting CC=gcc and CXX=g++

message("${CXX}")

#-std=c++0x to enable c++0x(c++11) suport, this is mandatory for Boost
#NDEBUG and BOOST_UBLAS_NDEBUG improve the Boost uBlas performance
#dramatically
set(CMAKE_CXX_FLAGS "-std=c++0x -DNDEBUG -DBOOST_UBLAS_NDEBUG")

#the multithread version of boost should be used. This should be 
#the default setting
find_package(Boost 1.51 COMPONENTS system filesystem regex graph timer serialization REQUIRED)

#mpi for parallel computing
find_package(MPI)
if (MPI_CXX_FOUND)
  find_package(Boost 1.51 COMPONENTS mpi REQUIRED)
  add_definitions( -DMPI_SUPPORT )
  include_directories(${MPI_CXX_INCLUDE_PATH})
else()
  message( "MPI is not found, no multiprocess parallel computing supported" )
endif()

message("MPI checking done")

#OpenMP is only optional
find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    add_definitions( -DOpenMP_SUPPORT )
else()
    message("OpenMP is not found, no multithreading computing supported")
endif()

include_directories(${Boost_INCLUDE_DIRS})

add_subdirectory(utils)
add_subdirectory(mcl)
add_subdirectory(lrw)


