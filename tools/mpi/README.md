# MPI Ports Tester

This tool tests the MPI Ports feature of the MPI implementation.
It exchanges the endpoints using the file `.address` in the working directory.

## To build

```
$ mkdir build
$ cd build
$ cmake ..
$ make
```

Optionally set the `CXX` environment variable to the desired MPI compiler wrapper prior to building.

## To test

```
$ mpirun -np 2 ./requester &
$ mpirun -np 2 ./acceptor
```
