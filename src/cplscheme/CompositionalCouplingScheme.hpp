// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_COMPOSITIONALCOUPLINGSCHEME_HPP_
#define PRECICE_CPLSCHEME_COMPOSITIONALCOUPLINGSCHEME_HPP_

#include "CouplingScheme.hpp"
#include "Constants.hpp"
#include "SharedPointer.hpp"

#include "com/Communication.hpp"
#include "mesh/SharedPointer.hpp"
#include "tarch/logging/Log.h"

#include <vector>
#include <list>

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace cplscheme {

/**
 * @brief Acts as one coupling scheme, but consists of several composed ones.
 *
 * The sequence of how the schemes are added by addCouplingScheme() determines
 * their execution sequence. The schemes are executed from first added to last
 * added, following the rule that implicit schemes need to be converged before
 * executing the next explicit scheme/s.
 * An example with implicit schemes I1, ..., I5 and explicit
 * schemes E1, ..., E4 is given here. The schemes are added in the sequence
 *
 *   I1, I2, E1, E2, I3, E3, I4, I5, E4
 *
 * The execution is as follows:
 *
 *   1. I1 and I2 iterate until convergence
 *   2. E1 and E2 are computed
 *   3. I3 iterates until convergence
 *   4. E3 is computed
 *   5. I4 and I5 iterate until convergence
 *   6. E4 is computed
 *
 * Steps 2 and 3 are computed directly when 1 is converged, i.e., no additional
 * solver iteration/computation is necessary. Step 3, however, will require
 * additional solver computations for every additional iteration. Similarly,
 * 4 and 5 are computed directly when 3 is converged and 6 is computed directly
 * when 5 is converged.
 *
 * The grouping into such "active" sets is repeatedly done, in order to know
 * which coupling schemes should be executed. The implicit coupling scheme
 * actions constants::actionWrite/ReadIterationCheckpoint() is used to find out
 * about whether a coupling scheme does perform iterations.
 *
 * When, in an active set of implicit coupling schemes, a scheme is converged but
 * others not, the converged scheme is put on hold until all schemes in the
 * active set are converged. This assumes that the converged state of a coupling
 * scheme is not deteriorated by the further iteration of the
 * remaining non-converged schemes.
 *
 * The execution of each scheme does not only depend on this sequence, but also
 * on how the participants are configured to be first and second in the schemes.
 * If not configured properly, a deadlock might be created.
 */
class CompositionalCouplingScheme : public CouplingScheme
{
public:

  /**
   * @brief Constructor.
   */
  CompositionalCouplingScheme();

  /**
   * @brief Destructor, empty.
   */
  virtual ~CompositionalCouplingScheme() {}

  /**
   * @brief Adds another coupling scheme in parallel to this scheme.
   *
   * If this coupling scheme is a normal coupling scheme, an object of
   * CompositionalCouplingScheme will be created that contains this and the
   * new scheme in parallel. If this coupling scheme is already a composed
   * scheme, the new scheme will be added as another parallel scheme.
   *
   * @return Pointer to composition of coupling schemes.
   */
  void addCouplingScheme(PtrCouplingScheme scheme);

  /**
   * @brief Initializes the coupling scheme and establishes a communiation
   *        connection to the coupling partner.
   */
  virtual void initialize (
    double startTime,
    int    startTimesteps );

  /**
   * @brief Returns true, if initialize has been called.
   */
  virtual bool isInitialized() const;

  /**
   * @brief Initializes the data for first implicit coupling scheme iteration.
   *
   * Has to be called after initialize() and before advance().
   */
  virtual void initializeData();

  /**
   * @brief Adds newly computed time. Has to be called before every advance.
   */
  virtual void addComputedTime(double timeToAdd);

  /**
   * @brief Exchanges data and updates the state of the coupling scheme.
   */
  virtual void advance();

  /**
   * @brief Finalizes the coupling and disconnects communication.
   */
  virtual void finalize();

  /*
   * @brief returns list of all coupling partners
   */
  virtual std::vector<std::string> getCouplingPartners() const;

  /**
   * @brief Returns true, if data will be exchanged when calling advance().
   *
   * Also returns true after the last call of advance() at the end of the
   * simulation.
   *
   * @param lastSolverTimestepLength [IN] The length of the last timestep
   *        computed by the solver calling willDataBeExchanged().
   */
  virtual bool willDataBeExchanged(double lastSolverTimestepLength) const;

  /**
   * @brief Returns true, if data has been exchanged in last call of advance().
   */
  virtual bool hasDataBeenExchanged() const;

  /**
   * @brief Returns the currently computed time of the coupling scheme.
   *
   * This time is the minimum time of any coupling scheme in the composition.
   */
  virtual double getTime() const;

  /**
   * @brief Returns the currently computed timesteps of the coupling scheme.
   *
   * The timestep is the minimum timestep in any coupling scheme in the composition.
   */
  virtual int getTimesteps() const;

  /**
   * @brief Returns the maximal time to be computed.
   *
   * This is the maximum max time of the coupling schemes in the composition.
   */
  virtual double getMaxTime() const;

  /**
   * @brief Returns the maximal timesteps to be computed.
   */
  virtual int getMaxTimesteps() const;

  /**
   * @brief Returns current subiteration number in timestep.
   */
  //virtual int getSubIteration() const;

  /**
   * @brief Returns true, if timestep length is prescribed by the cpl scheme.
   *
   * If any of the solvers in the composition has a timestep length limit, this
   * counts as limit.
   */
  virtual bool hasTimestepLength() const;

  /**
   * @brief Returns the timestep length, if one is given by the coupling scheme.
   *
   * An assertion is thrown, if no valid timestep is given. Check with
   * hasTimestepLength().
   *
   * The smallest timestep length limit in the coupling scheme composition has
   * to be obeyed.
   */
  virtual double getTimestepLength() const;

  /**
   * @brief Returns the remaining timestep length of the current time step.
   *
   * This is not necessarily the timestep length limit the solver has to obeye
   * which is returned by getNextTimestepMaxLength().
   *
   * If no timestep length is precribed by the coupling scheme, always 0.0 is
   * returned.
   *
   * The maximum remainer of all composed coupling schemes is returned.
   */
  virtual double getThisTimestepRemainder() const;

  /**
   * @brief Returns part of the current timestep that has been computed already.
   *
   * This is the minimum of all computed timestep parts of the composed coupling
   * schemes.
   */
  virtual double getComputedTimestepPart() const;

  /**
   * @brief Returns the maximal length of the next timestep to be computed.
   *
   * If no timestep length is prescribed by the coupling scheme, always the
   * maximal double accuracy floating point number value is returned.
   *
   * This is the minimum of all max lengths of the composed coupling schemes.
   */
  virtual double getNextTimestepMaxLength() const;

  /**
   * @brief Returns true, when the coupled simulation is still ongoing.
   *
   * As long as one composed coupling scheme is still ongoing, returns true.
   */
  virtual bool isCouplingOngoing() const;

  /**
   * @brief Returns true, when the accessor can advance to the next timestep.
   *
   * Only true, if all composed coupling schemes have completed the timestep.
   */
  virtual bool isCouplingTimestepComplete() const;

  /**
   * @brief Returns true, if the given action has to be performed by the accessor.
   *
   * True, if any of the composed coupling schemes requires the action.
   */
  virtual bool isActionRequired(const std::string& actionName) const;

  /**
   * @brief Tells the coupling scheme that the accessor has performed the given
   *        action.
   */
  virtual void performedAction(const std::string& actionName);

  /**
   * @brief Returns the checkpointing timestep interval.
   *
   * Returns the smallest interval of all composed schemes, although it makes
   * not much sense to have different intervals.
   */
  virtual int getCheckpointTimestepInterval() const;

  /**
   * @brief Sets an action required to be performed by the accessor.
   */
  virtual void requireAction(const std::string& actionName);

  /**
   * @brief Returns a string representation of the current coupling state.
   */
  virtual std::string printCouplingState() const;

  /**
   * @brief Exports the state of the coupling scheme to file/s.
   *
   * Used for checkpointing.
   */
  virtual void exportState(const std::string& filenamePrefix) const;

  /**
   * @brief Imports the state of the coupling scheme from file/s.
   *
   * Used for checkpointing.
   */
  virtual void importState(const std::string& filenamePrefix);

  /**
   * @brief Send the state of the coupling scheme to another remote scheme.
   *
   * Used in client-server approach for parallel solvers. There, the solver
   * interface does hold a coupling scheme with no data but state. The state
   * is transferred between the solver coupling scheme and the server coupling
   * scheme via sendState and receiveState.
   */
  virtual void sendState (
   com::Communication::SharedPointer communication,
   int                   rankReceiver );

  /**
   * @brief Receive the state of the coupling scheme from another remote scheme.
   *
   * Used in client-server approach for parallel solvers. There, the solver
   * interface does hold a coupling scheme with no data but state. The state
   * is transferred between the solver coupling scheme and the server coupling
   * scheme via sendState and receiveState.
   */
  virtual void receiveState (
   com::Communication::SharedPointer communication,
   int                   rankSender );

private:

  /**
   * @brief Groups a coupling scheme with additional associated variables.
   */
  struct Scheme {
    // @brief The actual coupling scheme
    PtrCouplingScheme scheme;

    // @brief Excludes converged implicit schemes from some operations.
    //
    // When several implicit schemes are iterating their point of convergence
    // will be different in general. Converged schemes are automatically advanced
    // to the next timestep and are put on hold until all schemes are converged.
    //
    // This assumes that a once converged scheme does not leave the convergence
    // region again.
    bool onHold;

    /**
     * @brief Constructor.
     */
    Scheme(PtrCouplingScheme scheme)
    : scheme(scheme), onHold(false) {}
  };

  typedef std::list<Scheme> Schemes;
  typedef std::list<Scheme>::iterator SchemesIt;
  //typedef std::list<PtrCouplingScheme>::const_iterator ConstSchemesIt;

  // @brief Logging device.
  static tarch::logging::Log _log;

  // @brief Coupling schemes to be executed in parallel.
  Schemes _couplingSchemes;

  //Schemes _activeCouplingSchemes;

  // @brief Iterator to begin of coupling schemes currently active.
  SchemesIt _activeSchemesBegin;

  // @brief Iterator to behind the end of coupling schemes currently active.
  SchemesIt _activeSchemesEnd;

  // @brief Stores time added since last call of advance.
  double _lastAddedTime;

  /**
   * @brief Determines the current set of active coupling schemes.
   *
   * Is called in advance, after all active sub-schemes are advanced.
   *
   * @return True, if active schemes have changed and should be treated within
   *         the same call of advance().
   */
  bool determineActiveCouplingSchemes();

  /**
   * @brief Advances the active set of coupling schemes.
   *
   * Helper function for determineActiveCouplingSchemes().
   */
  void advanceActiveCouplingSchemes();
};

}} // namespace precice, cplscheme

#endif /* PRECICE_CPLSCHEME_COMPOSITIONALCOUPLINGSCHEME_HPP_ */
