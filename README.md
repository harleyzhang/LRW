## Synopsis

This software implements Limited Random Walk (LRW) graph clustering algorithm described in 

H. Zhang, J. Raitoharju, S. Kiranyaz, and M. Gabbouj, “Limited random walk algorithm for big graph data clustering,” Journal of Big Data, vol. 3, no. 1, p. 26, 2016.

Graph clustering is an important task for knowledge minning. LRW algorithm lets agents randomly walk on a graph and collect the probability of each agent land on a node. LRW can solve graph clustering problem with high accuracy. It uses parellel paradigm and are capable of clustering graphs with millions of vertices and hundres of millions of edges.

If you find this software useful, please cite the aforementioned paper. 

## Installation

This software requires cmake and boost libraries. OpenMP and MPI must be installed if parallel computing is required.

This software has been tested on Linux environment. 

Build with the following scripts

```bash

mkdir build

cd build

cmake ..

make 

```


Note: 

** Note that cmake 2.8 has problem to link against boost libraries that are not installed on system directories. This bug is fixed in cmake 3.01

** To build the source code, C++11x must be supported. Use compile option
-std=c++0x

** Must use gnu compiler, CLang does not support C++11 very well. Use gcc version
4.9.0 or above


** use NMP_NUM_THREADS to change the number of thread that the process spouse. 

** set compiler: 
  currently, there are bug in cmake 2.8 that you can use SET(CMAKE_C_COMPILER mpicc), it will be fixed in cmake 3.1.xx for the moment, you have to do
  export CC=mpicc
  export CXX=mpicxx


## Tests

test folder contains simple bash script of using this software


## License

LGPL v3

