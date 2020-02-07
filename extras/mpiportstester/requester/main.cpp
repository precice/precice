// Preamble {{{
// =============================================================================
//       @file main.cpp
// -----------------------------------------------------------------------------
//     @author Alexander Shukaev <"alexander" "." "shukaev" "@" "tum" "." "de">
// -----------------------------------------------------------------------------
// @maintainer Alexander Shukaev <"alexander" "." "shukaev" "@" "tum" "." "de">
// -----------------------------------------------------------------------------
//  @copyright Copyright (C) 2015,
//             Alexander Shukaev <"alexander" "." "shukaev" "@" "tum" "." "de">.
//             All rights reserved.
// -----------------------------------------------------------------------------
//    @license This program is free software: you can redistribute it and/or
//             modify it under the terms of the GNU General Public License as
//             published by the Free Software Foundation, either version 3 of
//             the License, or (at your option) any later version.
//
//             This program is distributed in the hope that it will be useful,
//             but WITHOUT ANY WARRANTY; without even the implied warranty of
//             MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//             General Public License for more details.
//
//             You should have received a copy of the GNU General Public License
//             along with this program.  If not, see
//             <http://www.gnu.org/licenses/>.
// =============================================================================
// }}} Preamble

#include <mpi.h>

#include <fstream>
#include <iostream>

using std::cout;
using std::ifstream;
using std::remove;
using std::rename;
using std::string;

int size = 0;
int rank = -1;
int ack  = 0;

struct Sentinel {
  Sentinel()
  {
    ::ack = 1;
  }

  ~Sentinel()
  {
    if (::rank != 0)
      return;

    cout << "Requester: " << (::ack == 4 ? "Success!" : "Failure!") << '\n';
  }
} sentinel;

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  MPI_Comm_size(MPI_COMM_WORLD, &::size);
  MPI_Comm_rank(MPI_COMM_WORLD, &::rank);

  if (::rank == 0) {
    char port_name[MPI_MAX_PORT_NAME];

    string address_file_name(".address");

    {
      ifstream ifs;

      do {
        ifs.open(address_file_name);
      } while (not ifs);

      ifs.getline(port_name, MPI_MAX_PORT_NAME);
    }

    cout << "Address: " << port_name << '\n';

    MPI_Comm communicator;

    MPI_Comm_connect(port_name, MPI_INFO_NULL, 0, MPI_COMM_SELF, &communicator);

    MPI_Recv(&::ack, 1, MPI_INT, 0, 0x11, communicator, MPI_STATUS_IGNORE);

    ::ack *= 2;

    MPI_Send(&::ack, 1, MPI_INT, 0, 0x22, communicator);
    MPI_Recv(&::ack, 1, MPI_INT, 0, 0x33, communicator, MPI_STATUS_IGNORE);
  }

  MPI_Barrier(MPI_COMM_WORLD);

  MPI_Finalize();
}
