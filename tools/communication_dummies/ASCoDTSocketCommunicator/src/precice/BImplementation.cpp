#include <mpi.h>

#include "precice/BImplementation.h"

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
  static int vertexes[] = { 1, 0, 1, 1, 0, 2, 2, 4, 4, 4 };

  return vector<int>(vertexes,
                     vertexes + sizeof (vertexes) / sizeof (*vertexes));
}

vector<double> getData() {
  int rank = getRank();

  static double  data_0[] = { 0.0, 0.0 };
  static double  data_1[] = { 0.0, 0.0, 0.0 };
  static double  data_2[] = { 0.0, 0.0 };
  static double* data_3;
  static double  data_4[] = { 0.0, 0.0, 0.0 };

  static double* data[] = { data_0, data_1, data_2, data_3, data_4 };
  static int     size[] = { sizeof (data_0) / sizeof (*data_0),
                            sizeof (data_1) / sizeof (*data_1),
                            sizeof (data_2) / sizeof (*data_2),
                            0,
                            sizeof (data_4) / sizeof (*data_4) };

  return vector<double>(data[rank], data[rank] + size[rank]);
}

void*
task(void* b) {
  ((BImplementation*)b)->task();
}
}

BImplementation::
BImplementation()
  : _rank(getRank()),
    _counter(0),
    _vertexes(getVertexes()),
    _addresses(gatherAddresses()),
    _data(getData()),
    _initialized(false) {
  pthread_mutex_init(&_task_mutex, 0);

  pthread_mutex_init(&_initialize_mutex, 0);
  pthread_cond_init(&_initialize_cond, 0);

  pthread_mutex_init(&_receive_mutex, 0);
  pthread_cond_init(&_receive_cond, 0);
}

BImplementation::~BImplementation() {
  pthread_mutex_destroy(&_task_mutex);

  pthread_mutex_destroy(&_initialize_mutex);
  pthread_cond_destroy(&_initialize_cond);

  pthread_mutex_destroy(&_receive_mutex);
  pthread_cond_destroy(&_receive_cond);
}

void
BImplementation::
main() {
  pthread_t thread;

  pthread_create(&thread, 0, precice::task, this);
  pthread_detach(thread);
}

void
BImplementation::
task() {
  pthread_mutex_lock(&_task_mutex);

  if (!initialize()) {
    return;
  }

  receive();

  for (int i = 0; i < _data.size(); ++i) {
    _data[i] += _rank + 1;
  }

  cout << "B" << _rank << endl;

  for (int i = 0; i < _data.size(); ++i) {
    cout << _data[i] << endl;
  }

  cout << endl;

  pthread_mutex_unlock(&_task_mutex);
}

bool
BImplementation::
initialize() {
  if (_rank == 0) {
    cout << "B"
         << _rank
         << ":"
         << " "
         << "Handshaking with A...";

    int tag = cantorPairing('B', 'A');

    if (_a) {
      _a->acknowledge('B', tag);
    }

    if (tag == cantorPairing('A', 'B')) {
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
BImplementation::
acknowledge(int identifier, int& tag) {
  cout << "B"
       << _rank
       << ":"
       << " "
       << "Acknowledging A...";

  if (tag == cantorPairing(identifier, 'B') && _a) {
    tag = cantorPairing('B', identifier);

    _a->initializeParallel(&_addresses[0], _addresses.size(),
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
BImplementation::
initialize(string const* addresses,
           int           addresses_size,
           int const*    vertexes,
           int           vertexes_size) {
  pthread_mutex_lock(&_initialize_mutex);

  _a_addresses.assign(addresses, addresses + addresses_size);
  _a_vertexes.assign(vertexes,   vertexes  + vertexes_size);

  {
    _initialized = true;

    pthread_cond_signal(&_initialize_cond);
  }

  pthread_mutex_unlock(&_initialize_mutex);
}

void
BImplementation::
receive() {
  pthread_mutex_lock(&_receive_mutex);

  while (_counter != _data.size()) {
    pthread_cond_wait(&_receive_cond, &_receive_mutex);
  }

  pthread_mutex_unlock(&_receive_mutex);
}

void
BImplementation::
setData(double data, int index, int a_rank, int& tag) {
  if (tag == cantorPairing(a_rank, _rank)) {
    tag = cantorPairing(_rank, a_rank);

    pthread_mutex_lock(&_receive_mutex);

    if (_counter == _data.size()) {
      _counter = 0;
    }

    _data[index] = data;

    if (_counter < _data.size()) {
      _counter++;
    }

    if (_counter == _data.size()) {
      pthread_cond_signal(&_receive_cond);
    }

    pthread_mutex_unlock(&_receive_mutex);
  } else {
    cout << "B"
         << _rank
         << ":"
         << " "
         << "Error"
         << ":"
         << " "
         << "Tag mismatch from"
         << " "
         << "A"
         << a_rank
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
