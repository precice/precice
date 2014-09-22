#include <mpi.h>

#include "precice/AImplementation.h"
#include "precice/CommunicatorCxx2SocketPlainPort.h"

#include "cantorPairing.h"
#include "gatherAddresses.h"

using std::cout;
using std::endl;
using std::string;
using std::vector;

namespace precice {
namespace {
int
getRank() {
  int rank = 0;

  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  return rank;
}

vector<int> getVertexes() {
  static int vertexes[] = { 0, 0, 1, 0, 1, 1, 2, 0, 1, 2 };

  return vector<int>(vertexes,
                     vertexes + sizeof (vertexes) / sizeof (*vertexes));
}

vector<double> getData() {
  int rank = getRank();

  static double data_0[] = { 10.0, 20.0, 40.0, 80.0 };
  static double data_1[] = { 30.0, 50.0, 60.0, 90.0 };
  static double data_2[] = { 70.0, 100.0 };

  static double* data[] = { data_0, data_1, data_2 };
  static int     size[] = { sizeof (data_0) / sizeof (*data_0),
                            sizeof (data_1) / sizeof (*data_1),
                            sizeof (data_2) / sizeof (*data_2) };

  return vector<double>(data[rank], data[rank] + size[rank]);
}

void*
task(void* a) {
  ((AImplementation*)a)->task();
}
}

AImplementation::
AImplementation()
  : _rank(getRank()),
    _counter(0),
    _vertexes(getVertexes()),
    _addresses(gatherAddresses()),
    _data(getData()),
    _initialized(false) {
  pthread_mutex_init(&_task_mutex, 0);

  pthread_mutex_init(&_initialize_mutex, 0);
  pthread_cond_init(&_initialize_cond, 0);
}

AImplementation::~AImplementation() {
  pthread_mutex_destroy(&_task_mutex);

  pthread_mutex_destroy(&_initialize_mutex);
  pthread_cond_destroy(&_initialize_cond);
}

void
AImplementation::
main() {
  pthread_t thread;

  pthread_create(&thread, 0, precice::task, this);
  pthread_detach(thread);
}

void
AImplementation::
task() {
  pthread_mutex_lock(&_task_mutex);

  if (!initialize()) {
    return;
  }

  send();

  pthread_mutex_unlock(&_task_mutex);
}

bool
AImplementation::
initialize() {
  if (_rank == 0) {
    cout << "A"
         << _rank
         << ":"
         << " "
         << "Handshaking with B...";

    int tag = cantorPairing('A', 'B');

    if (_b) {
      _b->acknowledge('A', tag);
    }

    if (tag == cantorPairing('B', 'A')) {
      cout << " "
           << "Success."
           << endl;
    } else {
      cout << " "
           << "Failure!"
           << endl;

      return false;
    }
  }

  pthread_mutex_lock(&_initialize_mutex);

  while (!_initialized) {
    pthread_cond_wait(&_initialize_cond, &_initialize_mutex);
  }

  _initialized = false;

  pthread_mutex_unlock(&_initialize_mutex);

  return true;
}

void
AImplementation::
acknowledge(int identifier, int& tag) {
  cout << "A"
       << _rank
       << ":"
       << " "
       << "Acknowledging B...";

  if (tag == cantorPairing(identifier, 'A') && _b) {
    tag = cantorPairing('A', identifier);

    _b->initializeParallel(&_addresses[0], _addresses.size(),
                           &_vertexes[0],  _vertexes.size());

    cout << " "
         << "Success."
         << endl;
  } else {
    cout << " "
         << "Failure!"
         << endl;
  }
}

void
AImplementation::
initialize(string const* addresses,
           int           addresses_size,
           int const*    vertexes,
           int           vertexes_size) {
  pthread_mutex_lock(&_initialize_mutex);

  _b_addresses.assign(addresses, addresses + addresses_size);
  _b_vertexes.assign(vertexes,   vertexes  + vertexes_size);

  {
    _initialized = true;

    pthread_cond_signal(&_initialize_cond);
  }

  pthread_mutex_unlock(&_initialize_mutex);
}

void
AImplementation::
send() {
  vector<int> b_indexes(_b_addresses.size(), 0);

  for (int i = 0, j = 0; i < _vertexes.size(); ++i) {
    int a_rank =   _vertexes[i];
    int b_rank = _b_vertexes[i];

    if (a_rank == _rank) {
      send(_data[j++], b_indexes[b_rank], b_rank);
    }

    b_indexes[b_rank]++;
  }
}

void
AImplementation::
send(double data, int index, int b_rank) {
  string address = _b_addresses[b_rank];
  string host    = address.substr(0, address.find(":"));
  string port    = address.substr(host.length() + 1,
                                  address.length() - host.length() - 1);

  cout << host << ":" << port << endl;

  int buffer_size = atoi(getenv("PRECICE_A_BUFFER_SIZE"));

  // CommunicatorCxx2SocketPlainPort* communicator =
  // new CommunicatorCxx2SocketPlainPort(host.c_str(),
  // atoi(port.c_str()),
  // buffer_size);
  CommunicatorCxx2SocketPlainPort communicator(host.c_str(),
                                               atoi(port.c_str()),
                                               buffer_size);

  int tag = cantorPairing(_rank, b_rank);

  communicator.setData(data, index, _rank, tag);

  if (tag != cantorPairing(b_rank, _rank)) {
    cout << "A"
         << _rank
         << ":"
         << " "
         << "Error"
         << ":"
         << " "
         << "Tag mismatch from"
         << " "
         << "B"
         << b_rank
         << "!"
         << endl;
  }
}
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

  MPI_Init_thread(&argc, &argv, MPI_THREAD_SERIALIZED, &provided);
#endif

#ifdef _WIN32
  MAIN_LOOP(false);
#else
  main_loop_(false);
#endif
}
