// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IMPL_ABSTRACTDATAACTION_HPP_
#define PRECICE_IMPL_ABSTRACTDATAACTION_HPP_

#include "MeshContext.hpp"
#include "SharedPointer.hpp"
#include "utils/Helpers.hpp"
#include <vector>

namespace precice {
namespace impl {

/**
 * @brief Interface and common functionliaty for operations on coupling data.
 */
class AbstractDataAction
{
public:

   /**
    * @brief Defines the timing of application of the data action.
    */
   enum TimingConstants {
      ALWAYS_PRIOR,             // Everytime, before advancing cpl scheme
      ALWAYS_POST,              // Everytime, after advancing cpl scheme
      ON_EXCHANGE_PRIOR,        // On data exchange, before advancing cpl scheme
      ON_EXCHANGE_POST,         // On data exchange, after advancing cpl scheme
      ON_ADVANCE_TIMESTEP_POST  // On advancing to next dt, after adv. cpl scheme
   };

   /**
    * @brief Constructor.
    */
   AbstractDataAction ( TimingConstants timing )
   : _timing(timing)
   {}

   /**
    * @brief Destructor.
    */
   virtual ~AbstractDataAction() {}

   virtual void performAction (
     double dt,
     double computedPartFullDt,
     double fullDt ) =0;

   TimingConstants getTiming () const
   {
      return _timing;
   }

private:

   TimingConstants _timing;
};


}} // namespace precice, impl

#endif /* PRECICE_IMPL_ABSTRACTDATAACTION_HPP_ */
