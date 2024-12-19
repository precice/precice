#include <benchmark/benchmark.h>
#include <iostream>
#include <mpi.h>

int main(int argc, char **argv)
{
  // Setup MPI on one rank
  MPI_Init(nullptr, nullptr);
  int size;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  if (size > 1) {
    std::cerr << "ERROR: running on multiple ranks is not supported!\n";
    return 1;
  }

  ::benchmark::Initialize(&argc, argv);
  if (::benchmark::ReportUnrecognizedArguments(argc, argv))
    return 1;
  ::benchmark::RunSpecifiedBenchmarks();
  ::benchmark::Shutdown();

  MPI_Finalize();
  return 0;
}
