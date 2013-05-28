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

class CompositionalCouplingScheme : public CouplingScheme
{
public:

   /**
    * @brief Constructor.
    *
    * @param initialCouplingScheme [IN] First coupling scheme in the composition.
    *        Used to initialize the compositional coupling scheme.
    */
   CompositionalCouplingScheme(PtrCouplingScheme initialCouplingScheme);

   /**
    * @brief Destructor, empty.
    */
   virtual ~CompositionalCouplingScheme() {}

   /**
    * @brief Adds a coupling scheme to be executed in parallel with others.
    */
   void addParallelCouplingScheme( PtrCouplingScheme couplingScheme );

   /**
    * @brief Initializes the coupling scheme.
    *
    * Sets up communication to coupling partner, initializes coupling state.
    */
   virtual void initialize (
     double startTime,
     int    startTimestep );

   /**
    * @brief Initializes the data for first implicit coupling scheme iteration.
    *
    * Has to be called after initialize() and before advance().
    */
   virtual void initializeData();

   /**
    * @brief Adds newly computed time. Has to be called before every advance.
    */
   void addComputedTime ( double timeToAdd );

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
   virtual std::vector<std::string> getCouplingPartners (
     const std::string& accessorName ) const;

   virtual void sendState (
     com::PtrCommunication communication,
     int                   rankReceiver );

   virtual void receiveState (
     com::PtrCommunication communication,
     int                   rankSender );

   virtual std::string printCouplingState() const;

   virtual void exportState(io::TXTWriter& writer) const;

   virtual void importState(io::TXTReader& reader);

private:

   // @brief Logging device.
   static tarch::logging::Log _log;

   // @brief Coupling schemes to be executed in parallel.
   std::vector<PtrCouplingScheme> _couplingSchemes;
};

}} // namespace precice, cplscheme

#endif /* PRECICE_CPLSCHEME_COMPOSITIONALCOUPLINGSCHEME_HPP_ */
