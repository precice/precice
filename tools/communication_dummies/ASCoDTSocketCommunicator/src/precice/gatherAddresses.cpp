#include <mpi.h>

#include "gatherAddresses.h"

#include <iostream>

using std::cout;
using std::endl;
using std::string;
using std::vector;

string
retrieveSocketAddress();

namespace precice {
vector<string>
gatherAddresses() {
  string address = retrieveSocketAddress();

  int size = 1;
  int rank = 0;

  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  vector<string> addresses;

  addresses.resize(size);

  addresses[rank] = address;

  if (rank == 0) {
    for (int i = 1; i < size; ++i) {
      MPI_Status status;

      int length = 0;

      MPI_Recv(&length, 1, MPI_INT, i, 2000, MPI_COMM_WORLD, &status);

      string address(length, '\0');

      // cout << "Length: " << length << endl;

      MPI_Recv(&address[0], length, MPI_CHAR, i, 2001, MPI_COMM_WORLD, &status);

      // cout << "Rank: " << i << endl << "Address: " << address << endl <<
      // endl;

      addresses[i] = address;
    }
  } else {
    int length = address.length();

    MPI_Send(&length, 1, MPI_INT, 0, 2000, MPI_COMM_WORLD);
    MPI_Send(&address[0], length, MPI_CHAR, 0, 2001, MPI_COMM_WORLD);
  }

  return addresses;
}
}
