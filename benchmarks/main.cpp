#include <benchmark/benchmark.h>

#ifndef PRECICE_NO_MPI
#include <iostream>
#include <mpi.h>
#endif

int main(int argc, char **argv)
{
#ifndef PRECICE_NO_MPI
  // Setup MPI on one rank
  MPI_Init(nullptr, nullptr);
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size > 1) {
    std::cerr << "ERROR: running on multiple ranks is not supported!\n";
    return 1;
  }
#endif

  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  ::benchmark::RunSpecifiedBenchmarks();
  ::benchmark::Shutdown();

#ifndef PRECICE_NO_MPI
  MPI_Finalize();
#endif
  return 0;
}
