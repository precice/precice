#pragma once

#include "BaseCouplingScheme.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace cplscheme {

/**
 * @brief Coupling scheme for parallel coupling, i.e. simultaneous execution of two coupled participants
 *
 * For more information, look into Benjamin's thesis, Section 3.5. 
 * https://mediatum.ub.tum.de/doc/1320661/document.pdf
 */ 	
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

  virtual void initialize(double startTime, int startTimestep);

  virtual void initializeData();

  virtual void advance();


protected:
  /// merges send and receive data into one map (for parallel acceleration)
  virtual void mergeData();

  /// Returns all data (receive and send)
  DataMap& getAllData()
  {
    PRECICE_ASSERT(!doesFirstStep(), "Only the second participant should do the acceleration." );
    return _allData;
  }


private:
  logging::Logger _log{"cplscheme::ParallelCouplingScheme"};
  
  /// Map from data ID -> all data (receive and send) with that ID
  DataMap _allData;

  virtual void explicitAdvance();

  virtual void implicitAdvance();
};

}}
