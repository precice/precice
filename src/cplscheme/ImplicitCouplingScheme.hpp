// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_IMPLICITCOUPLINGSCHEME_HPP_
#define PRECICE_CPLSCHEME_IMPLICITCOUPLINGSCHEME_HPP_

#include "BaseCouplingScheme.hpp"
#include "SharedPointer.hpp"
#include "Constants.hpp"
#include "impl/SharedPointer.hpp"
#include "io/TXTTableWriter.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/PropertyContainer.hpp"
#include "mesh/Mesh.hpp"
#include "tarch/logging/Log.h"
#include "utils/Helpers.hpp"
#include "tarch/la/DynamicColumnMatrix.h"
#include "boost/tuple/tuple.hpp"

namespace precice {
  namespace cplscheme {
    namespace tests {
      class ImplicitCouplingSchemeTest;
    }
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace cplscheme {

/**
 * Abstract class that provides the basic functionalities for implicit coupling,
 * i.e. subiterating in every timestep to converge towards the strong solution.
 * The functionalities that differ for the classcial serial coupling and the
 * parallel coupling are implemented in the subclasses SerialImplicitCouplingScheme
 * and ParallelImplicitCouplingScheme.
 * Please look at ./coupling_steering.pdf for a brief sketch of the differences
 * between the serial and parallel implicit coupling.
 *
 * @brief Abstract coupling scheme with iterations per timestep to achieve strong solution.
 */
class ImplicitCouplingScheme : public BaseCouplingScheme
{
public:

  /**
   * @brief Constructor.
   *
   *
   * @param maxTime [IN] Simulation time limit, or UNDEFINED_TIME.
   * @param maxTimesteps [IN] Simulation timestep limit, or UNDEFINED_TIMESTEPS.
   * @param timestepLength [IN] Simulation timestep length.
   * @param firstParticipant [IN] Name of first participant in coupling.
   * @param secondParticipant [IN] Name of second participant in coupling.
   * @param localParticipant [IN] Name of participant using this coupling scheme.
   * @param communication [IN] Communication object for com. between participants.
   * @param maxIterations [IN] Maximal iterations per coupling timestep.
   * @param monitorIterations [IN] If true, a txt file monitoring iterations is
   *                          written.
   */
  ImplicitCouplingScheme (
    double                maxTime,
    int                   maxTimesteps,
    double                timestepLength,
    int                   validDigits,
    const std::string&    firstParticipant,
    const std::string&    secondParticipant,
    const std::string&    localParticipant,
    com::PtrCommunication communication,
    int                   maxIterations,
    constants::TimesteppingMethod dtMethod);

  /**
   * @brief Sets order of predictor of interface values for first participant.
   *
   * The first participant in the implicit coupling scheme has to take some
   * initial guess for the interface values computed by the second participant.
   * In order to improve this initial guess, an extrapolation from previous
   * timesteps can be performed.
   *
   * The standard predictor is of order zero, i.e., simply the converged values
   * of the last timestep are taken as initial guess for the coupling iterations.
   * Currently, an order 1 predictor is implement besides that.
   */
  void setExtrapolationOrder ( int order );

  /**
   * @brief Adds a measure to determine the convergence of coupling iterations.
   */
  void addConvergenceMeasure (
    int                         dataID,
    bool                        suffices,
    impl::PtrConvergenceMeasure measure );

  /**
   * @brief Set a coupling iteration post-processing technique.
   */
  void setIterationPostProcessing ( impl::PtrPostProcessing postProcessing );

  /*
   * @brief returns list of all coupling partners
   */
  virtual std::vector<std::string> getCouplingPartners () const;

  virtual void sendState (
    com::PtrCommunication communication,
    int                   rankReceiver );

  virtual void receiveState (
    com::PtrCommunication communication,
    int                   rankSender );

  virtual std::string printCouplingState() const;

  virtual void exportState(const std::string& filenamePrefix) const;

  virtual void importState(const std::string& filenamePrefix);

protected:

  void setIterationToPlot(int iterationToPlot){
    _iterationToPlot = iterationToPlot;
  }

  void setTimestepToPlot(int timestepToPlot){
    _timestepToPlot = timestepToPlot;
  }

  void setTimeToPlot(double timeToPlot){
    _timeToPlot = timeToPlot;
  }

  void setIterations(int iterations){
    _iterations = iterations;
  }

  int getIterations(){
    return _iterations;
  }

  int getTotalIterations(){
    return _totalIterations;
  }

  void increaseIterations(){
    _iterations++;
  }

  void increaseTotalIterations(){
    _totalIterations++;
  }

  void increaseIterationToPlot(){
    _iterationToPlot++;
  }

  int getMaxIterations(){
    return _maxIterations;
  }

  int getExtrapolationOrder(){
    return _extrapolationOrder;
  }



  void newConvergenceMeasurements();

  /**
   * @brief Updates internal state of coupling scheme for next timestep.
   */
  void timestepCompleted();

  /**
   * @brief Updates the convergence measurement of local send data.
   */
  bool measureConvergence();

  void extrapolateData(DataMap& data);


private:

  typedef tarch::la::DynamicColumnMatrix<double> DataMatrix;

  typedef tarch::la::DynamicVector<double> DataVector;

  // @brief Logging device.
  static tarch::logging::Log _log;


  // @brief Writes residuals to file.
//  io::TXTTableWriter _residualWriterL1;
//  io::TXT_communicationTableWriter _residualWriterL2;

  // @brief Writes value amplification to file.
//  io::TXTTableWriter _amplificationWriter;

  typedef boost::tuple<int,impl::PtrConvergenceMeasure> MeasureTuple;



  // @brief Limit of iterations during one timestep.
  int _maxIterations;

  // @brief Number of iteration in current timestep.
  int _iterationToPlot;
  int _timestepToPlot;
  double _timeToPlot;

  // @brief Number of iterations in current timestep.
  int _iterations;

  // @brief Number of total iterations performed.
  int _totalIterations;

//  void writeResidual (
//    const utils::DynVector& values,
//    const utils::DynVector& oldValues );

  friend class tests::ImplicitCouplingSchemeTest;
};

}} // namespace precice, cplscheme

#endif /* PRECICE_CPLSCHEME_IMPLICITCOUPLINGSCHEME_HPP_ */
