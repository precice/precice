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
bool                  MasterSlave::_isPrimaryRank = false;
bool                  MasterSlave::_isSecondaryRank  = false;
com::PtrCommunication MasterSlave::_communication;

logging::Logger MasterSlave::_log("utils::MasterSlave");

void MasterSlave::configure(Rank rank, int size)
{
  PRECICE_TRACE(rank, size);
  _rank = rank;
  _size = size;
  PRECICE_ASSERT(_rank != -1 && _size != -1);
  _isPrimaryRank = (rank == 0) && _size != 1;
  _isSecondaryRank  = (rank != 0);
  PRECICE_DEBUG("isSecondaryRank: {}, isPrimaryRank: {}", _isSecondaryRank, _isPrimaryRank);
}

Rank MasterSlave::getRank()
{
  return _rank;
}

int MasterSlave::getSize()
{
  return _size;
}

bool MasterSlave::isPrimary()
{
  return _isPrimaryRank;
}

bool MasterSlave::isSecondary()
{
  return _isSecondaryRank;
}

bool MasterSlave::isParallel()
{
  return _isPrimaryRank || _isSecondaryRank;
}

double MasterSlave::l2norm(const Eigen::VectorXd &vec)
{
  PRECICE_TRACE();

  if (not _isPrimaryRank && not _isSecondaryRank) { //old case
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
  /* old loop over all secondary ranks solution
  if(_isSecondaryRank){
    _communication->send(localSum2, 0);
    _communication->receive(globalSum2, 0);
  }
  if(_isPrimaryRank){
    globalSum2 += localSum2;
    for(Rank secondaryRank = 1; secondaryRank < _size; secondaryRank++){
      _communication->receive(localSum2, secondaryRank);
      globalSum2 += localSum2;
    }
    for(Rank secondaryRank = 1; secondaryRank < _size; secondaryRank++){
      _communication->send(globalSum2, secondaryRank);
    }
  }
  */
  return sqrt(globalSum2);
}

double MasterSlave::dot(const Eigen::VectorXd &vec1, const Eigen::VectorXd &vec2)
{
  PRECICE_TRACE();

  if (not _isPrimaryRank && not _isSecondaryRank) { //old case
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

  // old loop over all secondary ranks solution
  /*
  if(_isSecondaryRank){
    _communication->send(localSum, 0);
    _communication->receive(globalSum, 0);
  }
  if(_isPrimaryRank){
    globalSum += localSum;
    for(Rank secondaryRank = 1; secondaryRank < _size; secondaryRank++){
      _communication->receive(localSum, secondaryRank);
      globalSum += localSum;
    }
    for(Rank secondaryRank = 1; secondaryRank < _size; secondaryRank++){
      _communication->send(globalSum, secondaryRank);
    }
  }
  */
  return globalSum;
}

void MasterSlave::reset()
{
  PRECICE_TRACE();
  _isPrimaryRank = false;
  _isSecondaryRank  = false;
  _rank     = -1;
  _size     = -1;
}

void MasterSlave::reduceSum(precice::span<const double> sendData, precice::span<double> rcvData)
{
  PRECICE_TRACE();

  if (not _isPrimaryRank && not _isSecondaryRank) {
    std::copy(sendData.begin(), sendData.end(), rcvData.begin());
    return;
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());

  if (_isSecondaryRank) {
    // send local result to primary rank
    _communication->reduceSum(sendData, rcvData, 0);
  }

  if (_isPrimaryRank) {
    // receive local results from secondary ranks, apply SUM
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

  if (not _isPrimaryRank && not _isSecondaryRank) {
    rcvData = sendData;
    return;
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());

  if (_isSecondaryRank) {
    // send local result to primary rank
    _communication->reduceSum(sendData, rcvData, 0);
  }

  if (_isPrimaryRank) {
    // receive local results from secondary ranks, apply SUM
    _communication->reduceSum(sendData, rcvData);
  }
}

void MasterSlave::allreduceSum(precice::span<const double> sendData, precice::span<double> rcvData)
{
  PRECICE_TRACE();

  if (not _isPrimaryRank && not _isSecondaryRank) {
    std::copy(sendData.begin(), sendData.end(), rcvData.begin());
    return;
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());

  if (_isSecondaryRank) {
    // send local result to primary rank, receive reduced result from primary rank
    _communication->allreduceSum(sendData, rcvData, 0);
  }

  if (_isPrimaryRank) {
    // receive local results from secondary ranks, apply SUM, send reduced result to secondary ranks
    _communication->allreduceSum(sendData, rcvData);
  }
}

void MasterSlave::allreduceSum(double &sendData, double &rcvData)
{
  PRECICE_TRACE();

  if (not _isPrimaryRank && not _isSecondaryRank) {
    rcvData = sendData;
    return;
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());

  if (_isSecondaryRank) {
    // send local result to primary rank, receive reduced result from primary rank
    _communication->allreduceSum(sendData, rcvData, 0);
  }

  if (_isPrimaryRank) {
    // receive local results from secondary ranks, apply SUM, send reduced result to secondary ranks
    _communication->allreduceSum(sendData, rcvData);
  }
}

void MasterSlave::allreduceSum(int &sendData, int &rcvData)
{
  PRECICE_TRACE();

  if (not _isPrimaryRank && not _isSecondaryRank) {
    rcvData = sendData;
    return;
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());

  if (_isSecondaryRank) {
    // send local result to primary rank, receive reduced result from primary rank
    _communication->allreduceSum(sendData, rcvData, 0);
  }

  if (_isPrimaryRank) {
    // receive local results from secondary ranks, apply SUM, send reduced result to secondary ranks
    _communication->allreduceSum(sendData, rcvData);
  }
}

void MasterSlave::broadcast(precice::span<double> values)
{
  PRECICE_TRACE();

  if (not _isPrimaryRank && not _isSecondaryRank) {
    return;
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());

  if (_isPrimaryRank) {
    // Broadcast (send) value.
    _communication->broadcast(values);
  }

  if (_isSecondaryRank) {
    // Broadcast (receive) value.
    _communication->broadcast(values, 0);
  }
}

void MasterSlave::broadcast(bool &value)
{
  PRECICE_TRACE();

  if (not _isPrimaryRank && not _isSecondaryRank) {
    return;
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());

  if (_isPrimaryRank) {
    // Broadcast (send) value.
    _communication->broadcast(value);
  }

  if (_isSecondaryRank) {
    // Broadcast (receive) value.
    _communication->broadcast(value, 0);
  }
}

void MasterSlave::broadcast(double &value)
{
  PRECICE_TRACE();

  if (not _isPrimaryRank && not _isSecondaryRank) {
    return;
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());

  if (_isPrimaryRank) {
    // Broadcast (send) value.
    _communication->broadcast(value);
  }

  if (_isSecondaryRank) {
    // Broadcast (receive) value.
    _communication->broadcast(value, 0);
  }
}

} // namespace utils
} // namespace precice
