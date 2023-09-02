#include "com/Communication.hpp"
#ifndef PRECICE_NO_MPI

#include "acceleration/impl/ParallelMatrixOperations.hpp"
#include "utils/IntraComm.hpp"

namespace precice::acceleration::impl {

void ParallelMatrixOperations::initialize(const bool needCyclicComm)
{
  PRECICE_TRACE();

  if (needCyclicComm && utils::IntraComm::isParallel()) {
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

  const auto size = utils::IntraComm::getSize();
  const auto rank = utils::IntraComm::getRank();
  com::connectCircularComm("IQNIMVJCyclicComm", "", rank, size, *cyclicCommLeft, *cyclicCommRight);

  _cyclicCommLeft  = std::move(cyclicCommLeft);
  _cyclicCommRight = std::move(cyclicCommRight);
}

void ParallelMatrixOperations::closeCircularCommunication()
{
  // Enforce consistency
  PRECICE_ASSERT(static_cast<bool>(_cyclicCommLeft) == static_cast<bool>(_cyclicCommRight));
  PRECICE_ASSERT(_needCyclicComm);

  if ((utils::IntraComm::getRank() % 2) == 0) {
    _cyclicCommLeft->closeConnection();
    _cyclicCommRight->closeConnection();
  } else {
    _cyclicCommRight->closeConnection();
    _cyclicCommLeft->closeConnection();
  }

  _cyclicCommRight = nullptr;
  _cyclicCommLeft  = nullptr;
}

} // namespace precice::acceleration::impl

#endif
