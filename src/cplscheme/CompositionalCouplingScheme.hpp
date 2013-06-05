// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_COMPOSITIONALCOUPLINGSCHEME_HPP_
#define PRECICE_CPLSCHEME_COMPOSITIONALCOUPLINGSCHEME_HPP_

#include "CouplingScheme.hpp"
#include "Constants.hpp"
#include "SharedPointer.hpp"
#include "com/SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"
#include "tarch/logging/Log.h"

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace cplscheme {

/**
 * @brief Acts as one coupling scheme, but consists of several composed ones.
 *
 */
class CompositionalCouplingScheme : public CouplingScheme
{
public:

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
     com::PtrCommunication communication,
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
     com::PtrCommunication communication,
     int                   rankSender );

private:

   // @brief Logging device.
   static tarch::logging::Log _log;

   // @brief Coupling schemes to be executed in parallel.
   std::vector<PtrCouplingScheme> _couplingSchemes;
};

}} // namespace precice, cplscheme

#endif /* PRECICE_CPLSCHEME_COMPOSITIONALCOUPLINGSCHEME_HPP_ */
