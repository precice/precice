#include <mpi.h>

#include <sstream>

#include "precice/BImplementation.h"
#include "precice/ReceiverCxx2SocketPlainPort.h"

#include "cantorPairing.h"
#include "gatherAddresses.h"

using std::cout;
using std::endl;
using std::ostringstream;
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
    _vertexes(getVertexes()),
    _addresses(gatherAddresses()),
    _data(getData()),
    _counter(0),
    _initialized(true),
    _received(true) {
  pthread_mutex_init(&_task_mutex, 0);

  pthread_mutex_init(&_initialize_mutex, 0);
  pthread_cond_init(&_initialize_cond, 0);

  pthread_mutex_init(&_receive_mutex, 0);
  pthread_cond_init(&_receive_cond, 0);
}

BImplementation::
~BImplementation() {
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
  pthread_mutex_lock(&_task_mutex); {
    // This method is blocking. When it returns `true`, it guarantees that the
    // connection with A has been established and that `_a_addresses` (addresses
    // of A) and `_a_vertexes` (vertexes of A) are initialized.
    if (!initialize()) {
      return;
    }

    // This method is blocking. When it returns, it guarantees that `_data` (all
    // the required data) is received from A and ready to be processed further.
    receive();

    process();

    // This method is blocking. When it returns `true`, it guarantees that
    // `_data` (all the required data) has not only been sent to A, but has
    // already been received by A!
    if (!send()) {
      return;
    }
  } pthread_mutex_unlock(&_task_mutex);
}

bool
BImplementation::
initialize() {
  if (_rank == 0) {
    ostringstream oss;

    oss << "B"
        << _rank
        << ":"
        << " "
        << "Handshaking with A...";

    int tag = cantorPairing('B', 'A');

    if (_a) {
      _a->acknowledge('B', tag);
    }

    if (tag == cantorPairing('A', 'B')) {
      oss << " "
          << "Success.";
    } else {
      oss << " "
          << "Failure!";

      // TODO: Implement `interruptInitialization()` and call
      // `interruptInitializationParallel()` here to unblock all the MPI
      // processes awaiting the condition below.
    }

    cout << oss.str() << endl;
  }

  pthread_mutex_lock(&_initialize_mutex); {
    _initialized = false;

    pthread_cond_signal(&_initialize_cond);

    while (!_initialized) {
      pthread_cond_wait(&_initialize_cond, &_initialize_mutex);
    }
  } pthread_mutex_unlock(&_initialize_mutex);

  return true;
}

void
BImplementation::
acknowledge(int identifier, int& tag) {
  ostringstream oss;

  oss << "B"
      << _rank
      << ":"
      << " "
      << "Acknowledging A...";

  if (tag == cantorPairing(identifier, 'B') && _a) {
    tag = cantorPairing('B', identifier);

    _a->initializeParallel(&_addresses[0], _addresses.size(),
                           &_vertexes[0],  _vertexes.size());

    oss << " "
        << "Success.";
  } else {
    oss << " "
        << "Failure!";
  }

  cout << oss.str() << endl;
}

void
BImplementation::
initialize(string const* addresses,
           int           addresses_size,
           int const*    vertexes,
           int           vertexes_size) {
  ostringstream oss;

  oss << "B"
      << _rank
      << ":"
      << " "
      << "Initializing...";

  pthread_mutex_lock(&_initialize_mutex); {
    while (_initialized) {
      pthread_cond_wait(&_initialize_cond, &_initialize_mutex);
    }

    _a_addresses.assign(addresses, addresses + addresses_size);
    _a_vertexes.assign(vertexes,   vertexes  + vertexes_size);

    _initialized = true;

    pthread_cond_signal(&_initialize_cond);
  } pthread_mutex_unlock(&_initialize_mutex);

  oss << " "
      << "Success.";

  cout << oss.str() << endl;
}

bool
BImplementation::
send() {
  bool sent = true;

  ostringstream oss;

  oss << "B"
      << _rank
      << ":"
      << " "
      << "Sending to A...";

  vector<int> a_indexes(_a_addresses.size(), 0);

  for (int i = 0, j = 0; i < _vertexes.size(); ++i) {
    int a_rank = _a_vertexes[i];
    int b_rank =   _vertexes[i];

    if (b_rank == _rank) {
      // This method is blocking. When it returns `true`, it guarantees that
      // data has not only been sent to A, but has already been received by A!
      sent &= send(_data[j++], a_indexes[a_rank], a_rank);
    }

    a_indexes[a_rank]++;
  }

  if (sent) {
    oss << " "
        << "Success.";
  } else {
    oss << " "
        << "Failure!";
  }

  cout << oss.str() << endl;

  return sent;
}

bool
BImplementation::
send(double data, int index, int a_rank) {
  string address = _a_addresses[a_rank];
  string host    = address.substr(0, address.find(":"));
  string port    = address.substr(host.length() + 1,
                                  address.length() - host.length() - 1);

  int buffer_size = atoi(getenv("PRECICE_B_BUFFER_SIZE"));

  ReceiverCxx2SocketPlainPort receiver(host.c_str(),
                                       atoi(port.c_str()),
                                       buffer_size);

  int tag = cantorPairing(_rank, a_rank);

  receiver.receive(data, index, _rank, tag);

  return tag == cantorPairing(a_rank, _rank);
}

void
BImplementation::
receive() {
  ostringstream oss;

  oss << "B"
      << _rank
      << ":"
      << " "
      << "Receiving from A...";

  pthread_mutex_lock(&_receive_mutex); {
    _counter  = 0;
    _received = (_data.size() == 0);

    pthread_cond_broadcast(&_receive_cond);

    // Block until we receive all the required data.
    // NOTE: Unblocked by `receive(double, int, int, int&)`.
    while (!_received) {
      pthread_cond_wait(&_receive_cond, &_receive_mutex);
    }
  } pthread_mutex_unlock(&_receive_mutex);

  oss << " "
      << "Success.";

  cout << oss.str() << endl;
}

void
BImplementation::
receive(double data, int index, int a_rank, int& tag) {
  ostringstream oss;

  oss << "B"
      << _rank
      << ":"
      << " "
      << "Receiving from"
      << " "
      << "A"
      << a_rank
      << "...";

  if (tag == cantorPairing(a_rank, _rank)) {
    tag = cantorPairing(_rank, a_rank);

    pthread_mutex_lock(&_receive_mutex); {
      // Block receiving if we have already received all the required data.
      // Please, understand that as a consequence this will also block sending
      // in A.
      // NOTE: Unblocked by `receive()`.
      while (_received) {
        pthread_cond_wait(&_receive_cond, &_receive_mutex);
      }

      if (index >= 0 && index < _data.size()) {
        _data[index] = data;

        if (_counter < _data.size()) {
          _counter++;
        }

        if (_counter == _data.size()) {
          _received = true;

          pthread_cond_signal(&_receive_cond);
        }

        oss << " "
            << "Success.";
      } else {
        oss << " "
            << "Failure!";
      }
    } pthread_mutex_unlock(&_receive_mutex);
  } else {
    // TODO: Implement `interruptReceive()` and call
    // `interruptReceiveParallel()` here to unblock all the MPI processes
    // awaiting `receive()`.

    oss << " "
        << "Failure!";
  }

  // cout << oss.str() << endl;
}

void
BImplementation::
process() {
  for (int i = 0; i < _data.size(); ++i) {
    _data[i] += _rank + 1;
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
