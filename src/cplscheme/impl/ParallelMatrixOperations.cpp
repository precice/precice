#ifndef PRECICE_NO_MPI

#include "ParallelMatrixOperations.hpp"
#include "utils/Globals.hpp"
#include "utils/MasterSlave.hpp"
#include "com/MPIPortsCommunication.hpp"
#include "com/Communication.hpp"
#include "logging/Logger.hpp"
#include <Eigen/Dense>


namespace precice {
namespace cplscheme {
namespace impl {

logging::Logger ParallelMatrixOperations::
_log("precice::cplscheme::impl::ParallelMatrixOperations");


ParallelMatrixOperations::ParallelMatrixOperations() :
  _cyclicCommLeft(nullptr),
  _cyclicCommRight(nullptr),
  _needCycliclComm(true)
{}

void ParallelMatrixOperations::initialize(
  com::Communication::SharedPointer leftComm,
  com::Communication::SharedPointer rightComm,
  bool needCyclicComm)
{
  TRACE();

  _needCycliclComm = needCyclicComm;
  if(utils::MasterSlave::_masterMode ||utils::MasterSlave::_slaveMode){

    _cyclicCommLeft = leftComm;
    _cyclicCommRight = rightComm;

    if(_needCycliclComm){
      assertion(_cyclicCommLeft.get() != NULL); assertion(_cyclicCommLeft->isConnected());
      assertion(_cyclicCommRight.get() != NULL); assertion(_cyclicCommRight->isConnected());
    }
  }
}




}}} // namespace precice, cplscheme, impl

#endif
