#pragma once

#include "BaseCouplingScheme.hpp"
#include "tarch/logging/Log.h"


namespace precice {
namespace cplscheme {

class MultiCouplingScheme : public BaseCouplingScheme
{
public:
  MultiCouplingScheme (
    double                maxTime,
    int                   maxTimesteps,
    double                timestepLength,
    int                   validDigits,
    const std::string&    localParticipant,
    std::vector<m2n::M2N::SharedPointer> communications,
    constants::TimesteppingMethod dtMethod,
    int                   maxIterations = 1)
    ;

  /// @brief Logging device.
  static tarch::logging::Log _log;

  virtual void initialize(double startTime, int startTimestep);

  virtual void initializeData();

  virtual void advance();

  /// @brief Adds data to be sent on data exchange and possibly be modified during coupling iterations.
  void addDataToSend (
    mesh::PtrData data,
    mesh::PtrMesh mesh,
    bool          initialize,
    int           index);

  /// @brief Adds data to be received on data exchange.
  void addDataToReceive (
    mesh::PtrData data,
    mesh::PtrMesh mesh,
    bool          initialize,
    int           index);

protected:
  /// @brief merges send and receive data into one map (for parallel post-processing)
  virtual void mergeData();

private:
  void sendData();
  void receiveData();
  void setupConvergenceMeasures();
  CouplingData* getData ( int dataID );

  /// @brief Communication device to the other coupling participant.
  std::vector<m2n::M2N::SharedPointer> _communications;

  /// @brief Map from data ID -> all data (receive and send) with that ID
  DataMap _allData;

  std::vector<DataMap> _receiveDataVector;
  std::vector<DataMap> _sendDataVector;


};

}}
