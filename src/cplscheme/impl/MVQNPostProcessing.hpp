// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License

#ifndef PRECICE_NO_MPI

#ifndef PRECICE_CPLSCHEME_MVQNPOSTPROCESSING_HPP_
#define PRECICE_CPLSCHEME_MVQNPOSTPROCESSING_HPP_

#include "BaseQNPostProcessing.hpp"
#include "mesh/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "tarch/la/DynamicColumnMatrix.h"
#include "tarch/la/DynamicMatrix.h"
#include "tarch/la/DynamicVector.h"
#include "com/Communication.hpp"
#include "io/TXTWriter.hpp"
#include "ParallelMatrixOperations.hpp"
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
class MVQNPostProcessing : public BaseQNPostProcessing
{
public:

  /**
   * @brief Constructor.
   */
   MVQNPostProcessing (
      double initialRelaxation,
      bool   forceInitialRelaxation,
      int    maxIterationsUsed,
      int    timestepsReused,
      int 	 filter,
      double singularityLimit,
      std::vector<int>    dataIDs,
      PtrPreconditioner preconditioner,
      bool   alwaysBuildJacobian);

   /**
    * @brief Destructor, empty.
    */
   virtual ~MVQNPostProcessing();


   /**
    * @brief Initializes the post-processing.
    */
   virtual void initialize(DataMap& cplData);


   /**
    * @brief Marks a iteration sequence as converged.
    *
    * called by the iterationsConverged() method in the BaseQNPostProcessing class
    * handles the postprocessing sepcific action after the convergence of one iteration
    */
   virtual void specializedIterationsConverged(DataMap& cplData);
  
private:

   /// @brief: stores the approximation of the inverse Jacobian of the system at current time step.
   Eigen::MatrixXd _invJacobian;

   /// @brief: stores the approximation of the inverse Jacobian from the previous time step.
   Eigen::MatrixXd _oldInvJacobian;

   /// @brief: stores the sub result (W-J_prev*V) for the current iteration
   Eigen::MatrixXd _Wtil;

   /// @brief: Communication between neighboring slaves, backwards
   com::Communication::SharedPointer _cyclicCommLeft;

   /// @brief: Communication between neighboring slaves, forward
   com::Communication::SharedPointer _cyclicCommRight;

   /// @brief: encapsulates matrix-matrix and matrix-vector multiplications for serial and parallel execution
   ParallelMatrixOperations _parMatrixOps;

   /** @brief: if true, the less efficient method to compute the quasi-Newton update is used,
   *   that explicitly builds the Jacobian in each iteration. If set to false this is only done
   *   in the very last iteration and the update is computed based on MATVEC products.
   */
   bool _alwaysBuildJacobian;

   /** @brief: comptes the MVQN update using QR decomposition of V,
    *        furthermore it updates the inverse of the system jacobian
    */
   virtual void computeQNUpdate(DataMap& cplData, Eigen::VectorXd& xUpdate);
   
   /// @brief: updates the V, W matrices (as well as the matrices for the secondary data)
   virtual void updateDifferenceMatrices(DataMap & cplData);

   /// @brief: computes underrelaxation for the secondary data
   virtual void computeUnderrelaxationSecondaryData(DataMap& cplData);
   
   /** @brief: computes the quasi-Newton update vector based on the matrices V and W using a QR
    *  decomposition of V. The decomposition is not re-computed en-block in every iteration
    *  but updated so that the new added column in V is incorporated in the decomposition.
    *
    *  This method rebuilds the Jacobian matrix and the matrix W_til in each iteration
    *  which is not necessary and thus inefficient.
    */
   void computeNewtonUpdate(DataMap& cplData, Eigen::VectorXd& update);
   
   /** @brief: computes the quasi-Newton update vector based on the same numerics as above.
    *  However, it exploits the fact that the matrix W_til can be updated according to V and W
    *  via the formula W_til.col(j) = W.col(j) - J_inv * V.col(j).
    *  Then, pure matrix-vector products are sufficient to compute the update within one iteration, i.e.,
    *  (1) x1 := J_prev*(-res) (2) y := Z(-res) (3) xUp := W_til*y + x1
    *  The Jacobian matrix only needs to be set up in the very last iteration of one time step, i.e.
    *  in iterationsConverged.
    */
   void computeNewtonUpdateEfficient(DataMap& cplData, Eigen::VectorXd& update);

   /** @brief: computes the pseudo inverse of V multiplied with V^T, i.e., Z = (V^TV)^-1V^T via QR-dec
    */
   void pseudoInverse(Eigen::MatrixXd& pseudoInverse);

   /** @brief: computes a explicit representation of the Jacobian, i.e., n x n matrix
    */
   void buildJacobian();

   /** @brief: re-computes the matrix _Wtil = ( W - J_prev * V) instead of updating it according to V
    */
   void buildWtil();

   // @brief: Removes one iteration from V,W matrices and adapts _matrixCols.
   virtual void removeMatrixColumn(int columnIndex);

};

}}} // namespace precice, cplscheme, impl

#endif /* PRECICE_CPLSCHEME_MVQNPOSTPROCESSING_HPP_ */
#endif /* PRECICE_NO_MPI */
