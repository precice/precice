#include "BiCouplingScheme.hpp"
#include "m2n/M2N.hpp"
#include "m2n/SharedPointer.hpp"
#include "utils/Helpers.hpp"

namespace precice {
namespace cplscheme {

BiCouplingScheme::BiCouplingScheme(
    double                        maxTime,
    int                           maxTimeWindows,
    double                        timeWindowSize,
    int                           validDigits,
    const std::string &           firstParticipant,
    const std::string &           secondParticipant,
    const std::string &           localParticipant,
    m2n::PtrM2N                   m2n,
    int                           maxIterations,
    CouplingMode                  cplMode,
    constants::TimesteppingMethod dtMethod)
    : BaseCouplingScheme(maxTime, maxTimeWindows, timeWindowSize, validDigits, firstParticipant,
                         secondParticipant, localParticipant, maxIterations, cplMode, dtMethod),
      _m2n(m2n)
{
}

void BiCouplingScheme::addDataToSend(
    mesh::PtrData data,
    mesh::PtrMesh mesh,
    bool          requiresInitialization)
{
  PRECICE_TRACE();
  int id = data->getID();
  if (!utils::contained(id, _sendData)) {
    PtrCouplingData     ptrCplData(new CouplingData(&(data->values()), mesh, requiresInitialization, data->getDimensions()));
    DataMap::value_type pair = std::make_pair(id, ptrCplData);
    _sendData.insert(pair);
  } else {
    PRECICE_ERROR("Data \"" << data->getName() << "\" cannot be added twice for sending.");
  }
}

void BiCouplingScheme::addDataToReceive(
    mesh::PtrData data,
    mesh::PtrMesh mesh,
    bool          requiresInitialization)
{
  PRECICE_TRACE();
  int id = data->getID();
  if (!utils::contained(id, _receiveData)) {
    PtrCouplingData     ptrCplData(new CouplingData(&(data->values()), mesh, requiresInitialization, data->getDimensions()));
    DataMap::value_type pair = std::make_pair(id, ptrCplData);
    _receiveData.insert(pair);
  } else {
    PRECICE_ERROR("Data \"" << data->getName() << "\" cannot be added twice for receiving.");
  }
}

bool BiCouplingScheme::receiveConvergence()
{
  PRECICE_ASSERT(doesFirstStep(), "For convergence information the receiving participant is always the first one.");
  bool convergence;
  _m2n->receive(convergence);
  return convergence;
}

CouplingData *BiCouplingScheme::getSendData(
    int dataID)
{
  PRECICE_TRACE(dataID);
  DataMap::iterator iter = _sendData.find(dataID);
  if (iter != _sendData.end()) {
    return &(*(iter->second));
  }
  return nullptr;
}

CouplingData *BiCouplingScheme::getReceiveData(
    int dataID)
{
  PRECICE_TRACE(dataID);
  DataMap::iterator iter = _receiveData.find(dataID);
  if (iter != _receiveData.end()) {
    return &(*(iter->second));
  }
  return nullptr;
}

} // namespace cplscheme
} // namespace precice