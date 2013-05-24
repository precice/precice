// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_COUPLINGSCHEME_HPP_
#define PRECICE_CPLSCHEME_COUPLINGSCHEME_HPP_

#include "CouplingData.hpp"
#include "mesh/Data.hpp"
#include "com/SharedPointer.hpp"
#include "com/Constants.hpp"
#include "utils/PointerVector.hpp"
#include "tarch/logging/Log.h"
#include <set>
#include <map>
#include <string>
#include <ostream>
#include <vector>
#include <limits>

namespace precice {
  namespace io {
    class TXTWriter;
    class TXTReader;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace cplscheme {

/**
 * @brief Abstract base class for all coupling schemes.
 *
 * ! General description
 * A coupling scheme computes the actions to be done by the coupled participants
 * (solvers) in time. It provides interface functions to setup, advance and
 * shutdown the coupling scheme and interface functions to query the state of
 * the coupling scheme and required actions of the participants.
 *
 * ! Usage
 * -# create an object of a concrete coupling scheme class
 *    (ExplicitCouplingScheme, e.g.)
 * -# add all meshes holding data to the coupling scheme by addMesh()
 * -# configure the object by adding subclass specific information
 * -# start the coupling scheme with initialize(), where the name of the local
 *    participant, i.e. the participant using the coupling scheme object, is
 *    needed
 * -# retrieve necessary information about sent/received data and the state of
 *    the coupled simulation
 * -# query and fulfill required actions
 * -# compute data to be sent (possibly taking into account received data from
 *    initialize())
 * -# advance the coupling scheme with advance(); where the maximum timestep
 *    length needs to be obeyed
 * -# ....
 * -# when the method isCouplingOngoing() returns false, call finalize() to
 *    stop the coupling scheme
 */
class CouplingScheme
{
public:

  // @brief Does not define a time limit for the coupled simulation.
  static const double UNDEFINED_TIME;

  // @brief Does not define a timestep limit for the coupled simulation.
  static const int UNDEFINED_TIMESTEPS;

  // @brief To be used, when the coupling timestep length is determined
  //        dynamically during the coupling.
  static const double UNDEFINED_TIMESTEP_LENGTH;

   /**
    * @brief Constructor.
    */
   CouplingScheme (
     double maxTime,
     int    maxTimesteps,
     double timestepLength,
     int    validDigits );

   /**
    * @brief Destructor.
    */
   virtual ~CouplingScheme() {}

   /**
    * @brief Adds data to be sent on data exchange and possibly be modified during
    *        coupling iterations.
    */
   void addDataToSend (
      mesh::PtrData data,
      bool          initialize );

   /**
    * @brief Adds data to be received on data exchange.
    */
   void addDataToReceive (
      mesh::PtrData data,
      bool          initialize );

   /**
    * @brief Sets the checkpointing timestep interval.
    */
   void setCheckointTimestepInterval ( int timestepInterval );

   /**
    * @brief Initializes the coupling scheme and establishes a communiation
    *        connection to the coupling partner.
    */
   virtual void initialize (
      double startTime,
      int    startTimesteps ) =0;

   /**
    * @brief Returns true, if initialize has been called.
    */
   bool isInitialized() const;

   /**
    * @brief Initializes the data for first implicit coupling scheme iteration.
    *
    * Has to be called after initialize() and before advance().
    */
   virtual void initializeData() =0;

   /**
    * @brief Adds newly computed time. Has to be called before every advance.
    */
   virtual void addComputedTime ( double timeToAdd ) =0;

   /**
    * @brief Exchanges data and updates the state of the coupling scheme.
    */
   virtual void advance() =0;

   /**
    * @brief Finalizes the coupling and disconnects communication.
    */
   virtual void finalize() =0;

   /*
    * @brief checks whether data with specified id is used
   */
   bool isDataUsed ( int dataID );

   /*
    * @brief checks whether any data is used
   */
   bool isDataUsed();

   /*
    * @brief returns list of all coupling partners
   */
   virtual std::vector<std::string> getCouplingPartners (
     const std::string& accessorName ) const =0;

   /**
    * @brief Returns true, if data will be exchanged when calling advance().
    *
    * Also returns true after the last call of advance() at the end of the
    * simulation.
    *
    * @param lastSolverTimestepLength [IN] The length of the last timestep
    *        computed by the solver calling willDataBeExchanged().
    */
   bool willDataBeExchanged ( double lastSolverTimestepLength ) const;

   /**
    * @brief Returns true, if data has been exchanged in last call of advance().
    */
   bool hasDataBeenExchanged() const;

   /**
    * @brief Returns the currently computed time of the coupling scheme.
    */
   double getTime() const;

   /**
    * @brief Returns the currently computed timesteps of the coupling scheme.
    */
   int getTimesteps() const;

   /**
    * @brief Returns the maximal time to be computed.
    */
   double getMaxTime() const;

   /**
    * @brief Returns the maximal timesteps to be computed.
    */
   int getMaxTimesteps() const;

   /**
    * @brief Returns current subiteration number in timestep.
    */
   int getSubIteration() const
   {
     return _subIteration;
   }

   /**
    * @brief Returns true, if timestep length is prescribed by the cpl scheme.
    */
   bool hasTimestepLength() const;

   /**
    * @brief Returns the timestep length, if one is given by the coupling scheme.
    *
    * An assertion is thrown, if no valid timestep is given. Check with
    * hasTimestepLength().
    */
   double getTimestepLength() const;

   /**
    * @brief Returns the remaining timestep length of the current time step.
    *
    * If no timestep length is precribed by the coupling scheme, always 0.0 is
    * returned.
    */
   double getThisTimestepRemainder() const;

   /**
    * @brief Returns part of the current timestep that has been computed already.
    */
   double getComputedTimestepPart() const;

   /**
    * @brief Returns the maximal length of the next timestep to be computed.
    *
    * If no timestep length is prescribed by the coupling scheme, always the
    * maximal double accuracy floating point number value is returned.
    */
   double getNextTimestepMaxLength() const;

   /**
    * @brief Returns the number of valid digits when compare times.
    */
   int getValidDigits() const;

   /**
    * @brief Returns true, when the coupled simulation is still ongoing.
    */
   bool isCouplingOngoing() const;

   /**
    * @brief Returns true, when the accessor can advance to the next timestep.
    */
   bool isCouplingTimestepComplete() const;

   /**
    * @brief Returns true, if the given action has to be performed by the accessor.
    */
   bool isActionRequired ( const std::string& actionName ) const;

   /**
    * @brief Tells the coupling scheme that the accessor has performed the given
    *        action.
    */
   void performedAction ( const std::string& actionName );

   /**
    * @brief Returns the checkpointing timestep interval.
    */
   int getCheckpointTimestepInterval() const;

   /**
    * @brief Sets an action required to be performed by the accessor.
    */
   void requireAction ( const std::string& actionName );

   /**
    * @brief Returns a string representation of the current coupling state.
    */
   virtual std::string printCouplingState() const =0;

   virtual void exportState(io::TXTWriter& writer) const =0;

   virtual void importState(io::TXTReader& reader) = 0;

   virtual void sendState (
     com::PtrCommunication communication,
     int                   rankReceiver );

   virtual void receiveState (
     com::PtrCommunication communication,
     int                   rankSender );

protected:

   typedef std::map<int,CouplingData> DataMap;

   struct State {
     int id;
     std::string name;
   };

   /**
    * @brief Sends data sendDataIDs given in mapCouplingData with communication.
    */
   std::vector<int> sendData ( com::PtrCommunication communication );

   /**
    * @brief Receives data sendDataIDs given in mapCouplingData with communication.
    */
   std::vector<int> receiveData ( com::PtrCommunication communication );

   /**
    * @brief Returns all data to be sent.
    */
   const DataMap& getSendData() const
   {
      return _sendData;
   }

   const DataMap& getReceiveData() const
   {
      return _receiveData;
   }

   /**
    * @brief Returns all data to be sent.
    */
   DataMap& getSendData()
   {
      return _sendData;
   }

   DataMap& getReceiveData()
   {
      return _receiveData;
   }

   /**
    * @brief Sets the values
    */
   CouplingData* getSendData ( int dataID );

   /**
    * @brief Returns all data to be received with data ID as given.
    */
   CouplingData* getReceiveData ( int dataID );

   /**
    * @brief Sets value for computed timestep part.
    */
   void setComputedTimestepPart ( double computedTimestepPart );

   /**
    * @brief Sets flag to determine whether data has been exchanged in the last
    *        coupling iteration.
    */
   void setHasDataBeenExchanged ( bool hasDataBeenExchanged );

   /**
    * @brief Sets the compouted time of the coupling scheme.
    *
    * Used from subclasses and when a checkpoint has been read.
    */
   void setTime ( double time )
   {
      _time = time;
   }

   /**
    * @brief Sets the computed timesteps of the coupling scheme.
    *
    * Used from subclasses and when a checkpoint has been read.
    */
   void setTimesteps ( int timesteps )
   {
      _timesteps = timesteps;
   }

   void setTimestepLength ( double timestepLength )
   {
     _timestepLength = timestepLength;
   }

   void setSubIteration ( int subIteration )
   {
     _subIteration = subIteration;
   }

//   void setMaxLengthNextTimestep ( double limit )
//   {
//      _maxLengthNextTimestep = limit;
//   }

   void setIsCouplingTimestepComplete ( bool isCouplingTimestepComplete )
   {
      _isCouplingTimestepComplete = isCouplingTimestepComplete;
   }

   void setIsInitialized ( bool isInitialized )
   {
     _isInitialized = isInitialized;
   }

   /**
    * @brief If any required actions are open, an error message is issued.
    */
   void checkCompletenessRequiredActions();

   /**
    * @brief Returns a string representing the basic state w/o actions.
    */
   std::string printBasicState() const;

   /**
    * @brief As the version without parameters, but with changed timestep and time.
    *
    * This version is used by the ImplicitCouplingScheme at the moment, which
    * needs to use the last timestep in the plotting when the iterations of
    * a timestep are converged.
    */
   std::string printBasicState(
     int    timesteps,
     double time) const;

   /**
    * @brief Returns a string representing the required actions.
    */
   std::string printActionsState() const;

private:

   // @brief Logging device.
   static tarch::logging::Log _log;

   double _maxTime;

   int _maxTimesteps;

   double _timestepLength;

   int _validDigits;

   double _time;

   double _computedTimestepPart;

   int _timesteps;

   int _subIteration;

   int _checkpointTimestepInterval;

   bool _isCouplingOngoing;

   bool _isCouplingTimestepComplete;

   // @brief True, if data has been exchanged between solvers.
   bool _hasDataBeenExchanged;

   // @brief True, if coupling has been initialized.
   bool _isInitialized;

   std::set<std::string> _actions;

   // @brief Map from data ID -> all send data with that ID
   DataMap _sendData;

   // @brief Map from data ID -> all receive data with that ID
   DataMap _receiveData;
};

}} // namespace precice, cplscheme

#endif /* PRECICE_CPLSCHEME_COUPLINGSCHEME_HPP_ */
