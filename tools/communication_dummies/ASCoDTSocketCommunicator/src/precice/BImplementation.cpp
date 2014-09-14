#include <mpi.h>

#include "precice/BImplementation.h"

#include "gatherAddresses.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

precice::BImplementation::
BImplementation() {
  pthread_mutex_init(&_mutex, 0);

  _addresses = gatherAddresses();
}

precice::BImplementation::~BImplementation() {
  pthread_mutex_destroy(&_mutex);
}

extern "C" void
#ifdef _WIN32
MAIN_LOOP(bool joinable);
#else
main_loop_(bool joinable);
#endif

int
main(int argc, char** argv) {
#ifdef Parallel
  int provided;

  MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
#endif

#ifdef _WIN32
  MAIN_LOOP(false);
#else
  main_loop_(false);
#endif
}

void
precice::
BImplementation::
main() {
  if (_a) {
    _a->initializeAddressesParallel(&_addresses[0], _addresses.size());
  }
}

void
precice::
BImplementation::
initializeAddresses(string const* addresses,
                    int const     addresses_size) {
  _a_addresses.assign(addresses, addresses + addresses_size);

  // for (vector<string>::const_iterator i = _a_addresses.begin();
  // i != _a_addresses.end();
  // ++i) {
  // cout << *i << endl;
  // }
}

void
precice::
BImplementation::
initializeVertexes(int const* vertexes,
                   int const  vertexes_size) {}

void
precice::
BImplementation::
setData(double const* data, int const data_size) {
  pthread_mutex_lock(&_mutex);

  for (int i = 0; i < data_size; ++i) {
    cout << data[i] << endl;
  }

  int rank = 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  cout << "Rank: " << rank << endl;

  pthread_mutex_unlock(&_mutex);
}
