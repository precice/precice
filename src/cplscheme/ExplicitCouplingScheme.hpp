// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_EXPLICITCOUPLINGSCHEME_HPP_
#define PRECICE_CPLSCHEME_EXPLICITCOUPLINGSCHEME_HPP_

#include "BaseCouplingScheme.hpp"
#include "Constants.hpp"
#include "com/SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"
#include "tarch/logging/Log.h"

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace cplscheme {

class ExplicitCouplingScheme : public BaseCouplingScheme
{
public:

   /**
    * @brief Constructor.
    */
   ExplicitCouplingScheme (
     double                maxTime,
     int                   maxTimesteps,
     double                timestepLength,
     int                   validDigits,
     const std::string&    firstParticipant,
     const std::string&    secondParticipant,
     const std::string&    localParticipantName,
     com::PtrCommunication communication,
     constants::TimesteppingMethod dtMethod );

   /**
    * @brief Destructor, virtual because of virtual functions.
    */
   virtual ~ExplicitCouplingScheme();

   /**
    * @brief Initializes the coupling scheme.
    *
    * Sets up communication to coupling partner, initializes coupling state.
    */
   virtual void initialize (
     double startTime,
     int    startTimestep );

   /**
    * @brief Initializes data for first participant from second participant.
    *
    * If this method is not used, the first participant has zero initial values
    * for its read data, before receiving data in advance(). If non-zero values
    * are needed, this has to be configured in the coupling-scheme XML
    * exchange-data tags. A call of initializeData() is then mandatory for the
    * second participant. It has to be called after initialize() and before
    * advance(). The second participant has to write the initial data values
    * to preCICE after initialize() and before initializeData().
    */
   virtual void initializeData();

   /**
    * @brief Adds newly computed time. Has to be called before every advance.
    */
   //void addComputedTime ( double timeToAdd );

   /**
    * @brief Advances within the coupling scheme.
    */
   virtual void advance();

   /**
    * @brief Finalizes the coupling scheme.
    */
   virtual void finalize();

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

   virtual void exportState(const std::string& filenamePrefix) const {}

   virtual void importState(const std::string& filenamePrefix) {}

protected:

   /**
    * @return True, if local participant is the one starting the explicit scheme.
    */
   bool doesFirstStep(){
     return _doesFirstStep;
   }


   /**
    * @return Communication device to the other coupling participant.
    */
   com::PtrCommunication getCommunication(){
     return _communication;
   }


   void setHasToSendInitData(bool hasToSendInitData){
       _hasToSendInitData = hasToSendInitData;
     }

     void setHasToReceiveInitData(bool hasToReceiveInitData){
       _hasToReceiveInitData = hasToReceiveInitData;
     }

     bool hasToSendInitData(){
       return _hasToSendInitData;
     }

     bool hasToReceiveInitData(){
       return _hasToReceiveInitData;
     }

     bool participantReceivesDt(){
       return _participantReceivesDt;
     }

     bool participantSetsDt(){
       return _participantSetsDt;
     }

private:

   // @brief Logging device.
   static tarch::logging::Log _log;

   // @brief Name of participant starting the explicit coupling scheme.
   std::string _firstParticipant;

   // @brief Name of participant waiting for data from first participant.
   std::string _secondParticipant;

   // @brief True, if local participant is the one starting the explicit scheme.
   bool _doesFirstStep;

   // @brief Communication device to the other coupling participant.
   com::PtrCommunication _communication;

   bool _participantSetsDt;

   bool _participantReceivesDt;

   // @brief to carry initData information from initialize to initData
   bool _hasToSendInitData;

   // @brief to carry initData information from initialize to initData
   bool _hasToReceiveInitData;
};

}} // namespace precice, cplscheme

#endif /* PRECICE_CPLSCHEME_EXPLICITCOUPLINGSCHEME_HPP_ */
