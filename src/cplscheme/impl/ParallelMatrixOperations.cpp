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
  P_TRACE();

  _needCycliclComm = needCyclicComm;
  if (utils::MasterSlave::isMaster() || utils::MasterSlave::isSlave()) {

    _cyclicCommLeft  = leftComm;
    _cyclicCommRight = rightComm;

    if (_needCycliclComm) {
      P_ASSERT(_cyclicCommLeft.get() != NULL);
      P_ASSERT(_cyclicCommLeft->isConnected());
      P_ASSERT(_cyclicCommRight.get() != NULL);
      P_ASSERT(_cyclicCommRight->isConnected());
    }
  }
}
}
}
} // namespace precice, cplscheme, impl

#endif
