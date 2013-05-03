// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_ACTION_ACTION_HPP_
#define PRECICE_ACTION_ACTION_HPP_

#include "mesh/SharedPointer.hpp"

namespace precice {
namespace action {

/**
 * @brief Abstract base class for configurable actions on data and/or meshes.
 *
 * Actions are executed on call of precice::SolverInterface::initialize(),
 * ::initializeData(), and ::advance(). They can change meshes and in particular
 * data values.
 */
class Action
{
public:

   /**
    * @brief Defines the time and place of application of the action.
    */
   enum Timing {
      ALWAYS_PRIOR,             // Everytime, before advancing cpl scheme
      ALWAYS_POST,              // Everytime, after advancing cpl scheme
      ON_EXCHANGE_PRIOR,        // On data exchange, before advancing cpl scheme
      ON_EXCHANGE_POST,         // On data exchange, after advancing cpl scheme
      ON_TIMESTEP_COMPLETE_POST  // On advancing to next dt, after adv. cpl scheme
   };

   /**
    * @brief Constructor.
    */
   Action (
     Timing               timing,
     const mesh::PtrMesh& mesh )
   :
     _timing(timing),
     _mesh(mesh)
   {}

   /**
    * @brief Destructor, empty.
    */
   virtual ~Action() {}

   /**
    * @brief Performs the action, to be overwritten by subclasses.
    *
    * @param dt [IN] Length of last local timestep computed.
    * @param computedPartFullDt [IN] Sum of all local timesteps of current
    *        global timestep.
    * @param fullDt [IN] Current global timestep length.
    */
   virtual void performAction (
     double time,
     double dt,
     double computedPartFullDt,
     double fullDt ) =0;

   /**
    * @brief Returns the timing of the action.
    */
   Timing getTiming() const
   {
      return _timing;
   }

   /**
    * @brief Returns the mesh carrying the data used in the action.
    */
   const mesh::PtrMesh& getMesh() const
   {
     return _mesh;
   }

private:

   // @brief Determines when the action will be executed.
   Timing _timing;

   // @brief Mesh carrying the data used in the action.
   mesh::PtrMesh _mesh;
};


}} // namespace precice, action

#endif /* PRECICE_ACTION_ACTION_HPP_ */
