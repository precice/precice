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
  TRACE();

  _needCycliclComm = needCyclicComm;
  if (utils::MasterSlave::isMaster() || utils::MasterSlave::isSlave()) {

    _cyclicCommLeft  = leftComm;
    _cyclicCommRight = rightComm;

    if (_needCycliclComm) {
      assertion(_cyclicCommLeft.get() != NULL);
      assertion(_cyclicCommLeft->isConnected());
      assertion(_cyclicCommRight.get() != NULL);
      assertion(_cyclicCommRight->isConnected());
    }
  }
}
}
}
} // namespace precice, cplscheme, impl

#endif
