// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_BROYDENPOSTPROCESSING_HPP_
#define PRECICE_CPLSCHEME_BROYDENPOSTPROCESSING_HPP_

#include "BaseQNPostProcessing.hpp"
#include "mesh/SharedPointer.hpp"
#include "tarch/logging/Log.h"
#include "tarch/la/DynamicColumnMatrix.h"
#include "tarch/la/DynamicMatrix.h"
#include "tarch/la/DynamicVector.h"
#include "io/TXTWriter.hpp"
#include <deque>

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace cplscheme {
namespace impl {

/**
 * @brief Multi vector quasi-Newton update scheme 
 *
 * Performs a multi vector quasi-Newton to accelerate the convergence of implicit coupling
 * iterations. A multi Broyden update, together with the reuse of the approximate inverse 
 * Jacobian from the old time step are used to approximate the inverse Jacobian. After every
 * coupling iteration, the data values used are enhanced by the new coupling iterates.
 *
 * If more coupling data is present than used to compute the MVQN post-processing,
 * this data is relaxed using the same linear combination as computed for the
 * MVQN-related data. The data is called "secondary" henceforth and additional
 * old value and data matrices are needed for it.
 */
class BroydenPostProcessing : public BaseQNPostProcessing
{
public:

  /**
   * @brief Constructor.
   */
   BroydenPostProcessing (
      double initialRelaxation,
      bool forceInitialRelaxation,
      int    maxIterationsUsed,
      int    timestepsReused,
      int 	 filter,
      double singularityLimit,
      std::vector<int>    dataIDs,
      PtrPreconditioner preconditioner);

   /**
    * @brief Destructor, empty.
    */
   virtual ~BroydenPostProcessing() {}


   /**
    * @brief Initializes the post-processing.
    */
   virtual void initialize(DataMap& cplData);

   /**
    * @brief Performs one post-processing step.
    *
    * Has to be called after every implicit coupling iteration.
    */
   virtual void performPPSecondaryData(DataMap& cplData);

   /**
    * @brief Marks a iteration sequence as converged.
    *
    * called by the iterationsConverged() method in the BaseQNPostProcessing class
    * handles the postprocessing sepcific action after the convergence of one iteration
    */
   virtual void specializedIterationsConverged(DataMap& cplData);
  
private:

   // remove this ofter debugging, not useful
   // ---------------------------------------
   std::fstream f;
   //----------------------------------------
   
   // @brief stores the approximation of the inverse Jacobian of the system at current time step.
   Matrix _invJacobian;
   Matrix _oldInvJacobian;
   
   int _maxColumns;
   int _currentColumns;
   
  // @brief comptes the MVQN update using QR decomposition of V, 
  //        furthermore it updates the inverse of the system jacobian
   virtual void computeQNUpdate(DataMap& cplData, DataValues& xUpdate);
   
      // @brief updates the V, W matrices (as well as the matrices for the secondary data)
   virtual void updateDifferenceMatrices(DataMap & cplData);

   // @brief computes underrelaxation for the secondary data
   virtual void computeUnderrelaxationSecondaryData(DataMap& cplData);
   void computeNewtonFactorsQRDecomposition(DataMap& cplData, DataValues& update);
};

}}} // namespace precice, cplscheme, impl

#endif /* PRECICE_CPLSCHEME_BROYDENPOSTPROCESSING_HPP_ */

