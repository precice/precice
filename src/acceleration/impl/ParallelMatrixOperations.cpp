#ifndef PRECICE_NO_MPI

#include "acceleration/impl/ParallelMatrixOperations.hpp"
#include "utils/MasterSlave.hpp"

namespace precice {
namespace acceleration {
namespace impl {

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

} // namespace impl
} // namespace acceleration
} // namespace precice

#endif
