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
    * @brief Initializes the coupling scheme.
    *
    * Sets up communication to coupling partner, initializes coupling state.
    */
   virtual void initialize (
     double startTime,
     int    startTimestep ) = 0;

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
    * @return Communication device to the other coupling participant.
    */
   com::PtrCommunication getCommunication(){
     return _communication;
   }


private:

   // @brief Logging device.
   static tarch::logging::Log _log;

};

}} // namespace precice, cplscheme

#endif /* PRECICE_CPLSCHEME_EXPLICITCOUPLINGSCHEME_HPP_ */
