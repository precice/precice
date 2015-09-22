/*
 * MMPostProcessing.hpp
 *
 *  Created on: Sep 18, 2015
 *      Author: Klaudius Scheufele
 */

#ifndef MMPOSTPROCESSING_HPP_
#define MMPOSTPROCESSING_HPP_

#include "PostProcessing.hpp"
#include "SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"
#include "tarch/logging/Log.h"
#include "tarch/la/DynamicColumnMatrix.h"
#include "tarch/la/DynamicMatrix.h"
#include "tarch/la/DynamicVector.h"
#include "QRFactorization.hpp"
#include "Eigen/Dense"
#include <deque>
#include <fstream>
#include <string.h>


// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace cplscheme {
namespace impl {

/**
 * @brief Base Class for quasi-Newton post processing schemes
 *
 */
class MMPostProcessing : public PostProcessing
{
public:

  /**
   * @brief Constructor.
   */
	MMPostProcessing (
	  impl::PtrPostProcessing coarseModelOptimization,
      int    maxIterationsUsed,
      int    timestepsReused,
      int    filter,
      double singularityLimit,
      bool estimateJacobian,
      std::vector<int>    fineDataIDs,
      std::vector<int>    coarseDataIDs,
      std::map<int,double>    scalings);

   /**
    * @brief Destructor, empty.
    */
   virtual ~MMPostProcessing() {}

   /**
    * @brief Returns all MM involved fine model data IDs.
    */
   virtual std::vector<int> getFineDataIDs() const
   {
      return _fineDataIDs;
   }

   /**
    * @brief Returns all MM involved coarse model data IDs.
    */
   virtual std::vector<int> getCoarseDataIDs() const
   {
      return _coarseDataIDs;
   }


   /**
    * @brief Initializes the post-processing.
    */
   virtual void initialize(DataMap& cplData);

   /**
    * @brief Performs one post-processing step.
    *
    * Has to be called after every implicit coupling iteration.
    */
   virtual void performPostProcessing(DataMap& cplData);


   /**
    * @brief Marks a iteration sequence as converged.
    *
    * Since convergence measurements are done outside the post-processing, this
    * method has to be used to signalize convergence to the post-processing.
    */
   virtual void iterationsConverged(DataMap& cplData);


   /**
    * @brief sets the design specification we want to meet for the objective function,
    *     i. e., we want to solve for argmin_x ||R(x) - q||, with R(x) = H(x) - x
    *     Usually we want to solve for a fixed-point of H, thus solving for argmin_x ||R(x)||
    *     with q=0.
    */
   virtual void setDesignSpecification(Eigen::VectorXd& q);


   /**
    * @brief Sets whether the solver has to evaluate the coarse or the fine model representation
    * steers the coupling scheme and the post processing.
    */
   virtual void setNextModelToEvaluate(BaseCouplingScheme::ModelResolution& nextModel){
	   _nextModelToEvaluate = nextModel;
   }

   /**
    * @brief Exports the current state of the post-processing to a file.
    */
   virtual void exportState(io::TXTWriter& writer);

   /**
    * @brief Imports the last exported state of the post-processing from file.
    *
    * Is empty at the moment!!!
    */
   virtual void importState(io::TXTReader& reader);


   // delete this:
   virtual int getDeletedColumns();


protected:

   // @brief Logging device.
   static tarch::logging::Log _log;

   impl::PtrPostProcessing _coarseModelOptimization;

   // @brief Maximum number of old data iterations kept.
   int _maxIterationsUsed;

   // @brief Maximum number of old timesteps (with data values) kept.
   int _timestepsReused;

   // @brief Determines sensitivity when two matrix columns are considered equal.
   //
   // When during the QR decomposition of the V matrix a pivot element smaller
   // than the singularity limit is found, the matrix is considered to be singular
   // and the corresponding (older) iteration is removed.
   double _singularityLimit;

   /**
    * @brief sets the design specification we want to meet for the objective function,
    *     i. e., we want to solve for argmin_x ||R(x) - q||, with R(x) = H(x) - x
    *     Usually we want to solve for a fixed-point of H, thus solving for argmin_x ||R(x)||
    *     with q=0.
    */
   Eigen::VectorXd _designSpecification;
   Eigen::VectorXd _coarseModel_designSpecification;

   /**
    * @brief Sets whether the solver has to evaluate the coarse or the fine model representation
    * steers the coupling scheme and the post processing.
    */
   BaseCouplingScheme::ModelResolution _nextModelToEvaluate;

   /// @brief Data IDs of data to be involved in the MM algorithm.
   std::vector<int> _fineDataIDs;
   std::vector<int> _coarseDataIDs;
   std::vector<int> _dataIDs;

   /// @brief Data IDs of data not involved in MM coefficient computation.
   std::vector<int> _secondaryDataIDs;

   /** @brief Scales data by fixed value.
    *
    * When more than one data is used to compute the PP update, all data should have similar
    * magnitude, in order to be similarly important in the least-squares solution.
    */
   std::map<int,double> _scalings;

   // @brief Indicates the first iteration, where constant relaxation is used.
   bool _firstIteration;

   // @brief Indicates the first time step, where constant relaxation is used
   //        later, we replace the constant relaxation by a qN-update from last time step.
   bool _firstTimeStep;

   /** @brief Indicates whether the Jacobian is stored explicitly (multi-vector method) or
    *         if a matrix-vector update is used without explicit representation of the Jacobian.
    */
   bool _estimateJacobian;


   /// @brief current iteration residuals of fine model coupling data
   Eigen::VectorXd _fineResiduals;

   /// @brief current iteration residuals of coarse model coupling data
   Eigen::VectorXd _coarseResiduals;

   /// @brief difference between solver input and output of fine model from last time step
   Eigen::VectorXd _fineOldResiduals;

   /// @brief difference between solver input and output of coarse model from last time step
   Eigen::VectorXd _coarseOldResiduals;

   /// @brief Temporary used in performPostProcessing(). output fine model
   Eigen::VectorXd _outputFineModelScaled;

   /// @brief Temporary used in performPostProcessing().
   //Eigen::VectorXd _scaledOldValues;

   /// @brief Temporary used in performPostProcessing(), output coarse model
   Eigen::VectorXd _outputCoarseModelScaled;

   /// @brief Temporary used in performPostProcessing(), input for model evaluation
   Eigen::VectorXd _input_Xstar;

   /// @brief Temporary used in performPostProcessing().
   //Eigen::VectorXd _coarseScaledOldValues;

   /// @brief Stores residual deltas for the fine model response
   Eigen::MatrixXd _matrixF;

   /// @brief Stores residual deltas for the coarse model response
   Eigen::MatrixXd _matrixC;

   /// @brief The pseudo inverse of the manifold mapping matrix, only stored and updated
   ///        if _estimateJacobian is set to true.
   Eigen::MatrixXd _MMMappingMatrix;
   Eigen::MatrixXd _MMMappingMatrix_prev;


   /** @brief Indices (of columns in F, C matrices) of 1st iterations of timesteps.
    *
    * When old timesteps are reused (_timestepsReused > 0), the indices of the
    * first iteration of each timestep needs to be stored, such that, e.g., all
    * iterations of the last timestep, or one specific iteration that leads to
    * a singular matrix in the QR decomposition can be removed and tracked.
    */
   std::deque<int> _matrixCols;

   /** @brief only needed for the parallel master-slave mode. stores the local dimensions,
    *        i.e., the offsets in _invJacobian for all processors
    */
   std::vector<int> _dimOffsets;

   /// @ brief only debugging info, remove this:
   int its,tSteps;
   int deletedColumns;

   /** @brief filter method that is used to maintain good conditioning of the least-squares system
    *		Either of two types: QR1FILTER or QR2Filter
    */
   int _filter;

   /** @brief: computes number of cols in least squares system, i.e, number of cols in
    * 		  _matrixV, _matrixW, _qrV, etc..
    *		  This is necessary only for master-slave mode, when some procs do not have
    *		  any nodes on the coupling interface. In this case, the matrices are not
    * 		  constructed and we have no information about the number of cols. This info
    * 		  is needed for master-slave communication.
    * 		  Number of its =! _cols in general.
    */
   int getLSSystemCols();
   int getLSSystemRows();

   /// @brief updates the V, W matrices (as well as the matrices for the secondary data)
   void updateDifferenceMatrices(DataMap & cplData);

   /// @brief scales the data values with the predefined scaling factor
   void scaling(DataMap & cplData);

   /** @brief registers the new solution x_k+1 (x_star) from the coarse model optimization
    *         problem as new input data for the fine model evaluation step. This has to be done
    *         in each iteration that performs the coarse model optimization.
    */
   void registerSolutionCoarseModelOptimization(DataMap& cplData);

   /// reverts the scaling of the data values and overwrites the old values with the updated ones
   //void undoScaling(DataMap & cplData);

   /** @brief: computes/updates the design specification for the coarse model optimization problem
    *     	   i. e., q_k = c(x_k) - T_k * (f(x_k) - q), q = 0 is the fine model design specification
    */
   void computeCoarseModelDesignSpecifiaction();

   /// @brief computes the quasi-Newton update using the specified pp scheme (MVQN, IQNILS)
   void computeQNUpdate(DataMap& cplData, Eigen::VectorXd& xUpdate);

   /// @brief Removes one iteration from V,W matrices and adapts _matrixCols.
   void removeMatrixColumn(int columnIndex);

   // need to move that in a class/header that encapsulates the Eigen data types
   void shiftSetFirst(Eigen::MatrixXd& A, Eigen::VectorXd& v);
   void appendFront(Eigen::MatrixXd& A, Eigen::VectorXd& v);
   void MMPostProcessing::removeColumnFromMatrix(Eigen::MatrixXd& A, int col);

};

}}} // namespace precice, cplscheme, impl

#endif /* MMPOSTPROCESSING_HPP_ */
