#ifndef PRECICE_NO_MPI

#include "ParallelMatrixOperations.hpp"
#include "com/MPIPortsCommunication.hpp"
#include "utils/MasterSlave.hpp"

namespace precice
{
namespace cplscheme
{
namespace impl
{

void ParallelMatrixOperations::initialize(
    com::PtrCommunication leftComm,
    com::PtrCommunication rightComm,
    bool                  needCyclicComm)
{
  PRECICE_TRACE();

  _needCycliclComm = needCyclicComm;
  if (utils::MasterSlave::isMaster() || utils::MasterSlave::isSlave()) {

    _cyclicCommLeft  = leftComm;
    _cyclicCommRight = rightComm;

    if (_needCycliclComm) {
      PRECICE_ASSERT(_cyclicCommLeft.get() != NULL);
      PRECICE_ASSERT(_cyclicCommLeft->isConnected());
      PRECICE_ASSERT(_cyclicCommRight.get() != NULL);
      PRECICE_ASSERT(_cyclicCommRight->isConnected());
    }
  }
}
}
}
} // namespace precice, cplscheme, impl

#endif
