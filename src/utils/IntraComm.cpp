//#ifndef PRECICE_NO_MPI

#include <Eigen/Core>
#include <cmath>
#include <memory>
#include <ostream>
#include <string>

#include "IntraComm.hpp"
#include "com/Communication.hpp"
#include "logging/LogMacros.hpp"
#include "logging/Logger.hpp"
#include "precice/types.hpp"
#include "utils/assertion.hpp"
#include "utils/span_tools.hpp"

namespace precice {
namespace utils {

Rank                  IntraComm::_rank     = -1;
int                   IntraComm::_size     = -1;
bool                  IntraComm::_isPrimary = false;
bool                  IntraComm::_isSecondary  = false;
com::PtrCommunication IntraComm::_communication;

logging::Logger IntraComm::_log("utils::IntraComm");

void IntraComm::configure(Rank rank, int size)
{
  PRECICE_TRACE(rank, size);
  _rank = rank;
  _size = size;
  PRECICE_ASSERT(_rank != -1 && _size != -1);
  _isPrimary = (rank == 0) && _size != 1;
  _isSecondary  = (rank != 0);
  PRECICE_DEBUG("isSecondary: {}, isPrimary: {}", _isSecondary, _isPrimary);
}

Rank IntraComm::getRank()
{
  return _rank;
}

int IntraComm::getSize()
{
  return _size;
}

bool IntraComm::isPrimary()
{
  return _isPrimary;
}

bool IntraComm::isSecondary()
{
  return _isSecondary;
}

bool IntraComm::isParallel()
{
  return _isPrimary || _isSecondary;
}

double IntraComm::l2norm(const Eigen::VectorXd &vec)
{
  PRECICE_TRACE();

  if (not _isPrimary && not _isSecondary) { //old case
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
  if(_isSecondary){
    _communication->send(localSum2, 0);
    _communication->receive(globalSum2, 0);
  }
  if(_isPrimary){
    globalSum2 += localSum2;
    for(Rank rankSecondary = 1; rankSecondary < _size; rankSecondary++){
      _communication->receive(localSum2, rankSecondary);
      globalSum2 += localSum2;
    }
    for(Rank rankSecondary = 1; rankSecondary < _size; rankSecondary++){
      _communication->send(globalSum2, rankSecondary);
    }
  }
  */
  return sqrt(globalSum2);
}

double IntraComm::dot(const Eigen::VectorXd &vec1, const Eigen::VectorXd &vec2)
{
  PRECICE_TRACE();

  if (not _isPrimary && not _isSecondary) { //old case
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
  if(_isSecondary){
    _communication->send(localSum, 0);
    _communication->receive(globalSum, 0);
  }
  if(_isPrimary){
    globalSum += localSum;
    for(Rank rankSecondary = 1; rankSecondary < _size; rankSecondary++){
      _communication->receive(localSum, rankSecondary);
      globalSum += localSum;
    }
    for(Rank rankSecondary = 1; rankSecondary < _size; rankSecondary++){
      _communication->send(globalSum, rankSecondary);
    }
  }
  */
  return globalSum;
}

void IntraComm::reset()
{
  PRECICE_TRACE();
  _isPrimary = false;
  _isSecondary  = false;
  _rank     = -1;
  _size     = -1;
}

void IntraComm::reduceSum(precice::span<const double> sendData, precice::span<double> rcvData)
{
  PRECICE_TRACE();

  if (not _isPrimary && not _isSecondary) {
    std::copy(sendData.begin(), sendData.end(), rcvData.begin());
    return;
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());

  if (_isSecondary) {
    // send local result to master
    _communication->reduceSum(sendData, rcvData, 0);
  }

  if (_isPrimary) {
    // receive local results from slaves, apply SUM
    _communication->reduceSum(sendData, rcvData);
  }
}

void IntraComm::reduceSum(const double &sendData, double &rcvData)
{
  PRECICE_TRACE();
  reduceSum(precice::refToSpan<const double>(sendData),
            precice::refToSpan<double>(rcvData));
}

void IntraComm::reduceSum(const int &sendData, int &rcvData)
{
  PRECICE_TRACE();

  if (not _isPrimary && not _isSecondary) {
    rcvData = sendData;
    return;
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());

  if (_isSecondary) {
    // send local result to master
    _communication->reduceSum(sendData, rcvData, 0);
  }

  if (_isPrimary) {
    // receive local results from slaves, apply SUM
    _communication->reduceSum(sendData, rcvData);
  }
}

void IntraComm::allreduceSum(precice::span<const double> sendData, precice::span<double> rcvData)
{
  PRECICE_TRACE();

  if (not _isPrimary && not _isSecondary) {
    std::copy(sendData.begin(), sendData.end(), rcvData.begin());
    return;
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());

  if (_isSecondary) {
    // send local result to master, receive reduced result from master
    _communication->allreduceSum(sendData, rcvData, 0);
  }

  if (_isPrimary) {
    // receive local results from slaves, apply SUM, send reduced result to slaves
    _communication->allreduceSum(sendData, rcvData);
  }
}

void IntraComm::allreduceSum(double &sendData, double &rcvData)
{
  PRECICE_TRACE();

  if (not _isPrimary && not _isSecondary) {
    rcvData = sendData;
    return;
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());

  if (_isSecondary) {
    // send local result to master, receive reduced result from master
    _communication->allreduceSum(sendData, rcvData, 0);
  }

  if (_isPrimary) {
    // receive local results from slaves, apply SUM, send reduced result to slaves
    _communication->allreduceSum(sendData, rcvData);
  }
}

void IntraComm::allreduceSum(int &sendData, int &rcvData)
{
  PRECICE_TRACE();

  if (not _isPrimary && not _isSecondary) {
    rcvData = sendData;
    return;
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());

  if (_isSecondary) {
    // send local result to master, receive reduced result from master
    _communication->allreduceSum(sendData, rcvData, 0);
  }

  if (_isPrimary) {
    // receive local results from slaves, apply SUM, send reduced result to slaves
    _communication->allreduceSum(sendData, rcvData);
  }
}

void IntraComm::broadcast(precice::span<double> values)
{
  PRECICE_TRACE();

  if (not _isPrimary && not _isSecondary) {
    return;
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());

  if (_isPrimary) {
    // Broadcast (send) value.
    _communication->broadcast(values);
  }

  if (_isSecondary) {
    // Broadcast (receive) value.
    _communication->broadcast(values, 0);
  }
}

void IntraComm::broadcast(bool &value)
{
  PRECICE_TRACE();

  if (not _isPrimary && not _isSecondary) {
    return;
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());

  if (_isPrimary) {
    // Broadcast (send) value.
    _communication->broadcast(value);
  }

  if (_isSecondary) {
    // Broadcast (receive) value.
    _communication->broadcast(value, 0);
  }
}

void IntraComm::broadcast(double &value)
{
  PRECICE_TRACE();

  if (not _isPrimary && not _isSecondary) {
    return;
  }

  PRECICE_ASSERT(_communication.get() != nullptr);
  PRECICE_ASSERT(_communication->isConnected());

  if (_isPrimary) {
    // Broadcast (send) value.
    _communication->broadcast(value);
  }

  if (_isSecondary) {
    // Broadcast (receive) value.
    _communication->broadcast(value, 0);
  }
}

} // namespace utils
} // namespace precice
