#ifndef PRECICE_NO_MPI

#ifndef PRECICE_COM_MPI_REQUEST_HPP_
#define PRECICE_COM_MPI_REQUEST_HPP_

#include "Request.hpp"

#include <mpi.h>

namespace precice {
namespace com {
class MPIRequest : public Request {
public:
  MPIRequest(MPI_Request request);

  bool test();

  void wait();

private:
  MPI_Request _request;
};
}
} // namespace precice, com

#endif /* PRECICE_COM_MPI_REQUEST_HPP_ */

#endif // not PRECICE_NO_MPI
