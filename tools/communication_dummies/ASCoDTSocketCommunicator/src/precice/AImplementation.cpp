#include <mpi.h>

#include "precice/AImplementation.h"
#include "precice/CommunicatorCxx2SocketPlainPort.h"

#include "gatherAddresses.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

precice::AImplementation::
AImplementation() {
  pthread_mutex_init(&_mutex, 0);

  _addresses = gatherAddresses();
}

precice::AImplementation::~AImplementation() {
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
AImplementation::
main() {
  if (_b) {
    _b->initializeAddressesParallel(&_addresses[0], _addresses.size());
  }

  double data[] = { 10.0, 20.0, 40.0, 80.0 };

  send(_b_addresses[0], data, 4);
}

void
precice::
AImplementation::
initializeAddresses(string const* addresses,
                    int const     addresses_size) {
  _b_addresses.assign(addresses, addresses + addresses_size);

  // for (vector<string>::const_iterator i = _b_addresses.begin();
  // i != _b_addresses.end();
  // ++i) {
  // cout << *i << endl;
  // }
}

void
precice::
AImplementation::
initializeVertexes(int const* vertexes,
                   int const  vertexes_size) {}

void
precice::
AImplementation::
send(string const address, double const* data, int const data_size) {
  pthread_mutex_lock(&_mutex);

  string host = address.substr(0, address.find(":"));
  string port = address.substr(host.length() + 1,
                               address.length() - host.length() - 1);

  cout << host << ":" << port << endl;

  CommunicatorCxx2SocketPlainPort* communicator =
    new CommunicatorCxx2SocketPlainPort(const_cast<char*>(host.c_str()),
                                        atoi(port.c_str()),
                                        4096);

  communicator->setData(data, data_size);

  pthread_mutex_unlock(&_mutex);
}
