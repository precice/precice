#include <mpi.h>

#include <sstream>

#include "precice/AImplementation.h"
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

vector<double> getExpectedData() {
  int rank = getRank();

  static double data_0[] = { 12.0, 21.0, 42.0, 85.0 };
  static double data_1[] = { 32.0, 51.0, 63.0, 95.0 };
  static double data_2[] = { 73.0, 105.0 };

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

AImplementation::
~AImplementation() {
  pthread_mutex_destroy(&_task_mutex);

  pthread_mutex_destroy(&_initialize_mutex);
  pthread_cond_destroy(&_initialize_cond);

  pthread_mutex_destroy(&_receive_mutex);
  pthread_cond_destroy(&_receive_cond);
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
  pthread_mutex_lock(&_task_mutex); {
    // This method is blocking. When it returns `true`, it guarantees that the
    // connection with B has been established and that `_b_addresses` (addresses
    // of B) and `_b_vertexes` (vertexes of B) are initialized.
    if (!initialize()) {
      return;
    }

    _data = getData();

    // This method is blocking. When it returns `true`, it guarantees that
    // `_data` (all the required data) has not only been sent to B, but has
    // already been received by B!
    if (!send()) {
      return;
    }

    // This method is blocking. When it returns, it guarantees that `_data` (all
    // the required data) is received from B and ready to be processed further.
    receive();

    validate();
  } pthread_mutex_unlock(&_task_mutex);
}

bool
AImplementation::
initialize() {
  if (_rank == 0) {
    ostringstream oss;

    oss << "A"
        << _rank
        << ":"
        << " "
        << "Handshaking with B...";

    int tag = cantorPairing('A', 'B');

    if (_b) {
      _b->acknowledge('A', tag);
    }

    if (tag == cantorPairing('B', 'A')) {
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
AImplementation::
acknowledge(int identifier, int& tag) {
  ostringstream oss;

  oss << "A"
      << _rank
      << ":"
      << " "
      << "Acknowledging B...";

  if (tag == cantorPairing(identifier, 'A') && _b) {
    tag = cantorPairing('A', identifier);

    _b->initializeParallel(&_addresses[0], _addresses.size(),
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
AImplementation::
initialize(string const* addresses,
           int           addresses_size,
           int const*    vertexes,
           int           vertexes_size) {
  ostringstream oss;

  oss << "A"
      << _rank
      << ":"
      << " "
      << "Initializing...";

  pthread_mutex_lock(&_initialize_mutex); {
    while (_initialized) {
      pthread_cond_wait(&_initialize_cond, &_initialize_mutex);
    }

    _b_addresses.assign(addresses, addresses + addresses_size);
    _b_vertexes.assign(vertexes,   vertexes  + vertexes_size);

    _initialized = true;

    pthread_cond_signal(&_initialize_cond);
  } pthread_mutex_unlock(&_initialize_mutex);

  oss << " "
      << "Success.";

  cout << oss.str() << endl;
}

bool
AImplementation::
send() {
  bool sent = true;

  ostringstream oss;

  oss << "A"
      << _rank
      << ":"
      << " "
      << "Sending to B...";

  vector<int> b_indexes(_b_addresses.size(), 0);

  for (int i = 0, j = 0; i < _vertexes.size(); ++i) {
    int a_rank =   _vertexes[i];
    int b_rank = _b_vertexes[i];

    if (a_rank == _rank) {
      // This method is blocking. When it returns `true`, it guarantees that
      // data has not only been sent to B, but has already been received by B!
      sent &= send(_data[j++], b_indexes[b_rank], b_rank);
    }

    b_indexes[b_rank]++;
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
AImplementation::
send(double data, int index, int b_rank) {
  string address = _b_addresses[b_rank];
  string host    = address.substr(0, address.find(":"));
  string port    = address.substr(host.length() + 1,
                                  address.length() - host.length() - 1);

  int buffer_size = atoi(getenv("PRECICE_A_BUFFER_SIZE"));

  ReceiverCxx2SocketPlainPort receiver(host.c_str(),
                                       atoi(port.c_str()),
                                       buffer_size);

  int tag = cantorPairing(_rank, b_rank);

  receiver.receive(data, index, _rank, tag);

  return tag == cantorPairing(b_rank, _rank);
}

void
AImplementation::
receive() {
  ostringstream oss;

  oss << "A"
      << _rank
      << ":"
      << " "
      << "Receiving from B...";

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
AImplementation::
receive(double data, int index, int b_rank, int& tag) {
  ostringstream oss;

  oss << "A"
      << _rank
      << ":"
      << " "
      << "Receiving from"
      << " "
      << "B"
      << b_rank
      << "...";

  if (tag == cantorPairing(b_rank, _rank)) {
    tag = cantorPairing(_rank, b_rank);

    pthread_mutex_lock(&_receive_mutex); {
      // Block receiving if we have already received all the required data.
      // Please, understand that as a consequence this will also block sending
      // in B.
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

bool
AImplementation::
validate() {
  bool valid = true;

  ostringstream oss;

  oss << "A"
      << _rank
      << ":"
      << " "
      << "Validating...";

  vector<double> expected_data = getExpectedData();

  for (int i = 0; i < _data.size(); ++i) {
    valid &= (_data[i] == expected_data[i]);
  }

  if (valid) {
    oss << " "
        << "Success.";
  } else {
    oss << " "
        << "Failure!";
  }

  cout << oss.str() << endl;

  return valid;
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
