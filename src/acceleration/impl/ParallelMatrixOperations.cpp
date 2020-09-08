#ifndef PRECICE_NO_MPI

#include "acceleration/impl/ParallelMatrixOperations.hpp"
#include "utils/MasterSlave.hpp"

namespace precice {
namespace acceleration {
namespace impl {

void ParallelMatrixOperations::initialize(const bool needCyclicComm)
{
  PRECICE_TRACE();

  if (needCyclicComm && (utils::MasterSlave::isMaster() || utils::MasterSlave::isSlave())) {
    _needCyclicComm = true;
    establishCircularCommunication();
  } else {
    _needCyclicComm = false;
  }
}

ParallelMatrixOperations::~ParallelMatrixOperations()
{
  PRECICE_TRACE();

  if (_needCyclicComm) {
    closeCircularCommunication();
  }
}

void ParallelMatrixOperations::establishCircularCommunication()
{
  PRECICE_ASSERT(_needCyclicComm);
  PRECICE_ASSERT(_cyclicCommRight == nullptr);
  PRECICE_ASSERT(_cyclicCommLeft == nullptr);

  com::PtrCommunication cyclicCommLeft  = com::PtrCommunication(new com::MPIPortsCommunication("."));
  com::PtrCommunication cyclicCommRight = com::PtrCommunication(new com::MPIPortsCommunication("."));

  const auto size     = utils::MasterSlave::getSize();
  const auto rank     = utils::MasterSlave::getRank();
  const int  prevProc = (rank - 1 + size) % size;
  const int  nextProc = (rank + 1) % size;

  std::string prefix   = "MVQNCyclicComm";
  std::string prevName = prefix + std::to_string(prevProc);
  std::string thisName = prefix + std::to_string(rank);
  std::string nextName = prefix + std::to_string(nextProc);
  if ((rank % 2) == 0) {
    cyclicCommLeft->prepareEstablishment(prevName, thisName);
    cyclicCommLeft->acceptConnection(prevName, thisName, "", 0);
    cyclicCommLeft->cleanupEstablishment(prevName, thisName);

    cyclicCommRight->requestConnection(thisName, nextName, "", 0, 1);
  } else {
    cyclicCommRight->requestConnection(thisName, nextName, "", 0, 1);

    cyclicCommLeft->prepareEstablishment(prevName, thisName);
    cyclicCommLeft->acceptConnection(prevName, thisName, "", 0);
    cyclicCommLeft->cleanupEstablishment(prevName, thisName);
  }

  _cyclicCommLeft  = std::move(cyclicCommLeft);
  _cyclicCommRight = std::move(cyclicCommRight);
}

void ParallelMatrixOperations::closeCircularCommunication()
{
  // Enforce consistency
  PRECICE_ASSERT(static_cast<bool>(_cyclicCommLeft) == static_cast<bool>(_cyclicCommRight));
  PRECICE_ASSERT(_needCyclicComm);

  if ((utils::MasterSlave::getRank() % 2) == 0) {
    _cyclicCommLeft->closeConnection();
    _cyclicCommRight->closeConnection();
  } else {
    _cyclicCommRight->closeConnection();
    _cyclicCommLeft->closeConnection();
  }

  _cyclicCommRight = nullptr;
  _cyclicCommLeft  = nullptr;
}

} // namespace impl
} // namespace acceleration
} // namespace precice

#endif
