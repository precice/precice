#include <algorithm>
#include <map>
#include <memory>
#include <ostream>
#include <type_traits>
#include <utility>

#include "BiCouplingScheme.hpp"
#include "cplscheme/BaseCouplingScheme.hpp"
#include "cplscheme/CouplingData.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "logging/LogMacros.hpp"
#include "m2n/M2N.hpp"
#include "m2n/SharedPointer.hpp"
#include "mesh/Data.hpp"
#include "precice/types.hpp"
#include "utils/Helpers.hpp"

namespace precice::cplscheme {

BiCouplingScheme::BiCouplingScheme(
    double                        maxTime,
    int                           maxTimeWindows,
    double                        timeWindowSize,
    int                           validDigits,
    std::string                   firstParticipant,
    std::string                   secondParticipant,
    const std::string &           localParticipant,
    m2n::PtrM2N                   m2n,
    int                           maxIterations,
    CouplingMode                  cplMode,
    constants::TimesteppingMethod dtMethod,
    int                           extrapolationOrder)
    : BaseCouplingScheme(maxTime, maxTimeWindows, timeWindowSize, validDigits, localParticipant, maxIterations, cplMode, dtMethod, extrapolationOrder),
      _m2n(std::move(m2n)),
      _firstParticipant(std::move(firstParticipant)),
      _secondParticipant(std::move(secondParticipant))
{
  PRECICE_ASSERT(_firstParticipant != _secondParticipant,
                 "First participant and second participant must have different names.");
  if (localParticipant == _firstParticipant) {
    setDoesFirstStep(true);
  } else if (localParticipant == _secondParticipant) {
    setDoesFirstStep(false);
  } else {
    PRECICE_ERROR("Name of local participant \"{}\" does not match any participant specified for the coupling scheme.",
                  localParticipant);
  }
}

void BiCouplingScheme::addDataToSend(
    const mesh::PtrData &data,
    mesh::PtrMesh        mesh,
    bool                 requiresInitialization)
{
  PRECICE_TRACE();
  // @todo factor out into BaseCouplingScheme, function should create new CouplingData in _allData, if it does not exist and return the corresponding PtrCouplingData
  int             id = data->getID();
  PtrCouplingData aCplData;
  if (!utils::contained(id, _cplData)) { // data is not used by this coupling scheme yet, create new CouplingData
    if (isExplicitCouplingScheme()) {
      aCplData = std::make_shared<CouplingData>(data, std::move(mesh), requiresInitialization);
    } else {
      aCplData = std::make_shared<CouplingData>(data, std::move(mesh), requiresInitialization, getExtrapolationOrder());
    }
    _cplData.emplace(id, aCplData);
  } else { // data is already used by another exchange of this coupling scheme, use existing CouplingData
    aCplData = _cplData[id];
  }

  if (!utils::contained(id, _sendData)) {
    PRECICE_ASSERT(_sendData.count(id) == 0, "Key already exists!");
    _sendData.emplace(id, aCplData);
  } else {
    PRECICE_ERROR("Data \"{0}\" cannot be added twice for sending. Please remove any duplicate <exchange data=\"{0}\" .../> tags", data->getName());
  }
}

void BiCouplingScheme::addDataToReceive(
    const mesh::PtrData &data,
    mesh::PtrMesh        mesh,
    bool                 requiresInitialization)
{
  PRECICE_TRACE();
  // @todo factor out into BaseCouplingScheme, function should create new CouplingData in _allData, if it does not exist and return the corresponding PtrCouplingData
  int             id = data->getID();
  PtrCouplingData aCplData;
  if (!utils::contained(id, _cplData)) { // data is not used by this coupling scheme yet, create new CouplingData
    if (isExplicitCouplingScheme()) {
      aCplData = std::make_shared<CouplingData>(data, std::move(mesh), requiresInitialization);
    } else {
      aCplData = std::make_shared<CouplingData>(data, std::move(mesh), requiresInitialization, getExtrapolationOrder());
    }
    _cplData.emplace(id, aCplData);
  } else { // data is already used by another exchange of this coupling scheme, use existing CouplingData
    aCplData = _cplData[id];
  }

  if (!utils::contained(id, _receiveData)) {
    PRECICE_ASSERT(_receiveData.count(id) == 0, "Key already exists!");
    _receiveData.emplace(id, aCplData);
  } else {
    PRECICE_ERROR("Data \"{0}\" cannot be added twice for receiving. Please remove any duplicate <exchange data=\"{0}\" ... /> tags", data->getName());
  }
}

void BiCouplingScheme::determineInitialDataExchange()
{
  determineInitialSend(_sendData);
  determineInitialReceive(_receiveData);
}

std::vector<std::string> BiCouplingScheme::getCouplingPartners() const
{
  std::vector<std::string> partnerNames;
  // Add non-local participant
  if (doesFirstStep()) {
    partnerNames.push_back(_secondParticipant);
  } else {
    partnerNames.push_back(_firstParticipant);
  }
  return partnerNames;
}

DataMap &BiCouplingScheme::getSendData()
{
  return _sendData;
}

DataMap &BiCouplingScheme::getReceiveData()
{
  return _receiveData;
}

CouplingData *BiCouplingScheme::getSendData(
    DataID dataID)
{
  PRECICE_TRACE(dataID);
  DataMap::iterator iter = _sendData.find(dataID);
  if (iter != _sendData.end()) {
    return &(*(iter->second));
  }
  return nullptr;
}

CouplingData *BiCouplingScheme::getReceiveData(
    DataID dataID)
{
  PRECICE_TRACE(dataID);
  DataMap::iterator iter = _receiveData.find(dataID);
  if (iter != _receiveData.end()) {
    return &(*(iter->second));
  }
  return nullptr;
}

m2n::PtrM2N BiCouplingScheme::getM2N() const
{
  PRECICE_ASSERT(_m2n);
  return _m2n;
}

bool BiCouplingScheme::hasAnySendData()
{
  return not getSendData().empty();
}

bool BiCouplingScheme::hasSendData(DataID dataID)
{
  return getSendData(dataID) != nullptr;
}

void BiCouplingScheme::overwriteReceiveData(std::string dataName, double relativeDt)
{
  bool mustOverwrite = true;
  PRECICE_ASSERT(math::greaterEquals(relativeDt, time::Storage::WINDOW_START), relativeDt);
  PRECICE_ASSERT(math::greaterEquals(time::Storage::WINDOW_END, relativeDt), relativeDt);
  for (auto &receiveData : _receiveData | boost::adaptors::map_values) {
    if (receiveData->getDataName() == dataName) {
      receiveData->storeValuesAtTime(relativeDt, receiveData->values(), mustOverwrite);
      return;
    }
  }
  PRECICE_ASSERT(false, "Data with name not found", dataName);
}

void BiCouplingScheme::loadReceiveDataFromStorage(std::string dataName, double relativeDt)
{
  PRECICE_ASSERT(math::greaterEquals(relativeDt, time::Storage::WINDOW_START), relativeDt);
  PRECICE_ASSERT(math::greaterEquals(time::Storage::WINDOW_END, relativeDt), relativeDt);
  for (auto &receiveData : _receiveData | boost::adaptors::map_values) {
    if (receiveData->getDataName() == dataName) {
      receiveData->values() = receiveData->getValuesAtTime(relativeDt);
      return;
    }
  }
  PRECICE_ASSERT(false, "Data with name not found", dataName);
}

void BiCouplingScheme::storeSendValuesAtTime(double relativeDt)
{
  for (auto &data : getSendData() | boost::adaptors::map_values) {
    data->storeValuesAtTime(relativeDt, data->values());
  }
}

void BiCouplingScheme::initializeSendDataStorage()
{
  for (auto &data : getSendData() | boost::adaptors::map_values) {
    data->storeValuesAtTime(time::Storage::WINDOW_START, data->values());
    data->storeValuesAtTime(time::Storage::WINDOW_END, data->values());
  }
}

std::vector<double> BiCouplingScheme::getReceiveTimes(std::string dataName)
{
  auto times = std::vector<double>();
  for (auto &data : _receiveData | boost::adaptors::map_values) {
    if (data->getDataName() == dataName) {
      auto timesVec = data->getStoredTimesAscending();
      PRECICE_ASSERT(timesVec.size() > 0, timesVec.size());
      for (int i = 0; i < timesVec.size(); i++) {
        times.push_back(timesVec(i));
      }
      return times;
    }
  }
  PRECICE_DEBUG("No data with dataName {} found in receive data. Returning empty.", dataName);
  // PRECICE_ASSERT(false);  // Reasonable assertion, but the test Integration/Serial/WatchIntegralScaleAndNoScale actually uses read-data that is not receiving any data and this assertion gets triggered. See also https://github.com/precice/precice/pull/1526
  return times;
}

} // namespace precice::cplscheme
