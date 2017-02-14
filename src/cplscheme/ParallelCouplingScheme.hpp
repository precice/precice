#pragma once

#include "BaseCouplingScheme.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace cplscheme {

class ParallelCouplingScheme : public BaseCouplingScheme
{
public:
  ParallelCouplingScheme (
    double                maxTime,
    int                   maxTimesteps,
    double                timestepLength,
    int                   validDigits,
    const std::string&    firstParticipant,
    const std::string&    secondParticipant,
    const std::string&    localParticipant,
    m2n::PtrM2N           m2n,
    constants::TimesteppingMethod dtMethod,
    CouplingMode          cplMode,
    int                   maxIterations = 1);

  static logging::Logger _log;

  virtual void initialize(double startTime, int startTimestep);

  virtual void initializeData();

  virtual void advance();


protected:
  /// merges send and receive data into one map (for parallel post-processing)
  virtual void mergeData();

  /// Returns all data (receive and send)
  DataMap& getAllData()
    {
      assertion(!doesFirstStep(), "Only the second participant should do the post processing." );
      return _allData;
    }


private:
  /// @brief Map from data ID -> all data (receive and send) with that ID
  DataMap _allData;

  virtual void explicitAdvance();

  virtual void implicitAdvance();
};

}}
