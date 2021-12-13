//#ifndef PRECICE_NO_MPI

#include <Eigen/Core>
#include <cmath>
#include <memory>
#include <ostream>
#include <string>

#include "MasterSlave.hpp"
#include "com/Communication.hpp"
#include "logging/LogMacros.hpp"
#include "logging/Logger.hpp"
#include "precice/types.hpp"
#include "utils/assertion.hpp"
#include "utils/span_tools.hpp"

namespace precice {
namespace utils {

Rank                  MasterSlave::_rank     = -1;
int                   MasterSlave::_size     = -1;
bool                  MasterSlave::_isMaster = false;
bool                  MasterSlave::_isSlave  = false;
com::PtrCommunication MasterSlave::_communication;

logging::Logger MasterSlave::_log("utils::MasterSlave");

void MasterSlave::configure(Rank rank, int size)
{
  PRECICE_TRACE(rank, size);
  _rank = rank;
  _size = size;
  PRECICE_ASSERT(_rank != -1 && _size != -1);
  _isMaster = (rank == 0) && _size != 1;
  _isSlave  = (rank != 0);
  PRECICE_DEBUG("isSlave: {}, isMaster: {}", _isSlave, _isMaster);
}

Rank MasterSlave::getRank()
{
  return _rank;
}

int MasterSlave::getSize()
{
  return _size;
}

bool MasterSlave::isMaster()
{
  return _isMaster;
}

bool MasterSlave::isSlave()
{
  return _isSlave;
}

bool MasterSlave::isParallel()
{
  return _isMaster || _isSlave;
}

double MasterSlave::l2norm(const Eigen::VectorXd &vec)
{
  PRECICE_TRACE();

  if (not _isMaster && not _isSlave) { //old case
    return vec.norm();
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());
  double localSum2  = 0.0;
  double globalSum2 = 0.0;

  for (int i = 0; i < vec.size(); i++) {
    localSum2 += vec(i) * vec(i);
  }

  // localSum is modified, do not use afterwards
  allreduceSum(localSum2, globalSum2);
  /* old loop over all slaves solution
  if(_isSlave){
    _communication->send(localSum2, 0);
    _communication->receive(globalSum2, 0);
  }
  if(_isMaster){
    globalSum2 += localSum2;
    for(Rank rankSlave = 1; rankSlave < _size; rankSlave++){
      _communication->receive(localSum2, rankSlave);
      globalSum2 += localSum2;
    }
    for(Rank rankSlave = 1; rankSlave < _size; rankSlave++){
      _communication->send(globalSum2, rankSlave);
    }
  }
  */
  return sqrt(globalSum2);
}

double MasterSlave::dot(const Eigen::VectorXd &vec1, const Eigen::VectorXd &vec2)
{
  PRECICE_TRACE();

  if (not _isMaster && not _isSlave) { //old case
    return vec1.dot(vec2);
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());
  PRECICE_ASSERT(vec1.size() == vec2.size(), vec1.size(), vec2.size());
  double localSum  = 0.0;
  double globalSum = 0.0;

  for (int i = 0; i < vec1.size(); i++) {
    localSum += vec1(i) * vec2(i);
  }

  // localSum is modified, do not use afterwards
  allreduceSum(localSum, globalSum);

  // old loop over all slaves solution
  /*
  if(_isSlave){
    _communication->send(localSum, 0);
    _communication->receive(globalSum, 0);
  }
  if(_isMaster){
    globalSum += localSum;
    for(Rank rankSlave = 1; rankSlave < _size; rankSlave++){
      _communication->receive(localSum, rankSlave);
      globalSum += localSum;
    }
    for(Rank rankSlave = 1; rankSlave < _size; rankSlave++){
      _communication->send(globalSum, rankSlave);
    }
  }
  */
  return globalSum;
}

void MasterSlave::reset()
{
  PRECICE_TRACE();
  _isMaster = false;
  _isSlave  = false;
  _rank     = -1;
  _size     = -1;
}

void MasterSlave::reduceSum(precice::span<const double> sendData, precice::span<double> rcvData)
{
  PRECICE_TRACE();

  if (not _isMaster && not _isSlave) {
    std::copy(sendData.begin(), sendData.end(), rcvData.begin());
    return;
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());

  if (_isSlave) {
    // send local result to master
    _communication->reduceSum(sendData, rcvData, 0);
  }

  if (_isMaster) {
    // receive local results from slaves, apply SUM
    _communication->reduceSum(sendData, rcvData);
  }
}

void MasterSlave::reduceSum(const double &sendData, double &rcvData)
{
  PRECICE_TRACE();
  reduceSum(precice::refToSpan<const double>(sendData),
            precice::refToSpan<double>(rcvData));
}

void MasterSlave::reduceSum(const int &sendData, int &rcvData)
{
  PRECICE_TRACE();

  if (not _isMaster && not _isSlave) {
    rcvData = sendData;
    return;
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());

  if (_isSlave) {
    // send local result to master
    _communication->reduceSum(sendData, rcvData, 0);
  }

  if (_isMaster) {
    // receive local results from slaves, apply SUM
    _communication->reduceSum(sendData, rcvData);
  }
}

void MasterSlave::allreduceSum(precice::span<const double> sendData, precice::span<double> rcvData)
{
  PRECICE_TRACE();

  if (not _isMaster && not _isSlave) {
    std::copy(sendData.begin(), sendData.end(), rcvData.begin());
    return;
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());

  if (_isSlave) {
    // send local result to master, receive reduced result from master
    _communication->allreduceSum(sendData, rcvData, 0);
  }

  if (_isMaster) {
    // receive local results from slaves, apply SUM, send reduced result to slaves
    _communication->allreduceSum(sendData, rcvData);
  }
}

void MasterSlave::allreduceSum(double &sendData, double &rcvData)
{
  PRECICE_TRACE();

  if (not _isMaster && not _isSlave) {
    rcvData = sendData;
    return;
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());

  if (_isSlave) {
    // send local result to master, receive reduced result from master
    _communication->allreduceSum(sendData, rcvData, 0);
  }

  if (_isMaster) {
    // receive local results from slaves, apply SUM, send reduced result to slaves
    _communication->allreduceSum(sendData, rcvData);
  }
}

void MasterSlave::allreduceSum(int &sendData, int &rcvData)
{
  PRECICE_TRACE();

  if (not _isMaster && not _isSlave) {
    rcvData = sendData;
    return;
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());

  if (_isSlave) {
    // send local result to master, receive reduced result from master
    _communication->allreduceSum(sendData, rcvData, 0);
  }

  if (_isMaster) {
    // receive local results from slaves, apply SUM, send reduced result to slaves
    _communication->allreduceSum(sendData, rcvData);
  }
}

void MasterSlave::broadcast(precice::span<double> values)
{
  PRECICE_TRACE();

  if (not _isMaster && not _isSlave) {
    return;
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());

  if (_isMaster) {
    // Broadcast (send) value.
    _communication->broadcast(values);
  }

  if (_isSlave) {
    // Broadcast (receive) value.
    _communication->broadcast(values, 0);
  }
}

void MasterSlave::broadcast(bool &value)
{
  PRECICE_TRACE();

  if (not _isMaster && not _isSlave) {
    return;
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());

  if (_isMaster) {
    // Broadcast (send) value.
    _communication->broadcast(value);
  }

  if (_isSlave) {
    // Broadcast (receive) value.
    _communication->broadcast(value, 0);
  }
}

void MasterSlave::broadcast(double &value)
{
  PRECICE_TRACE();

  if (not _isMaster && not _isSlave) {
    return;
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());

  if (_isMaster) {
    // Broadcast (send) value.
    _communication->broadcast(value);
  }

  if (_isSlave) {
    // Broadcast (receive) value.
    _communication->broadcast(value, 0);
  }
}

} // namespace utils
} // namespace precice
