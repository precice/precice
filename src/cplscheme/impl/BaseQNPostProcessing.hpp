// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_BASEQNPOSTPROCESSING_HPP_
#define PRECICE_CPLSCHEME_BASEQNPOSTPROCESSING_HPP_

#include "PostProcessing.hpp"
#include "mesh/SharedPointer.hpp"
#include "tarch/logging/Log.h"
#include "tarch/la/DynamicColumnMatrix.h"
#include "tarch/la/DynamicMatrix.h"
#include "tarch/la/DynamicVector.h"
#include "QRFactorization.hpp"
#include "Preconditioner.hpp"
#include "Eigen/Dense"
#include <deque>
#include <fstream>
#include <string.h>



/* ****************************************************************************
 * 
 * A few comments concerning the present design choice.
 * 
 * All the functions from the base class BAseQNPostProcessing are specialized in
 * the sub classes as needed. This is done vi overwriting the base functions in 
 * the specialized sub classes and calling the respective base function after 
 * performing the specialized stuff (in order to perform the common, generalized 
 * computations i.e. handling of V,W matrices etc.)
 * However, for the performPostProcessing Method we decided (for better readability)
 * to have this method only in the base class, while introducing a function 
 * performPPSecondaryData that handles all the specialized stuff concerning post
 * processing for the secondary data in the sub classes.
 *
 * Another possibility would have been to introduce a bunch of functions like
 * initializeSpecialized(), removeMatrixColumnSpecialized(), 
 * iterationsConvergedSpecialized(), etc 
 * and call those function from the base class top down to the sub classes.
 *
 * The third possibility was to separate the approximation of the Jacobian from
 * the common stuff like handling V,W matrices in the post processing. 
 * Here, we have a class QNPostProcessing that handles the V,W stuff an d the basic
 * scheme of the QN update. Furthermore we have a base class (or rather interface) 
 * JacobianApproximation with sub classes MVQNAPX and IQNAPX that handle all the 
 * specialized stuff like Jacobian approximation, handling of secondary data etc. 
 * However, this approach is not feasible, as we have to call the function 
 * removeMatrixColumn() down in the specialized sub classes MVQNApx and IQNApx.
 * This is not possible as the function works on the V, W matrices that are 
 * completely treated by QNPostProcessing.
 * 
 * ****************************************************************************
 */




// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace cplscheme {
namespace impl {

/**
 * @brief Base Class for quasi-Newton post processing schemes
 * 
 */
class BaseQNPostProcessing : public PostProcessing
{
public:

  /**
   * @brief Constructor.
   */
   BaseQNPostProcessing (
      double initialRelaxation,
      bool   forceInitialRelaxation,
      int    maxIterationsUsed,
      int    timestepsReused,
      int    filter,
      double singularityLimit,
      std::vector<int>    dataIDs,
      PtrPreconditioner preconditioner);

   /**
    * @brief Destructor, empty.
    */
   virtual ~BaseQNPostProcessing() {
     
      //if ( _timingStream.is_open() ) {
      //  _timingStream.close ();
     // }
  }

   /**
    * @brief Returns all IQN involved data IDs.
    */
   virtual std::vector<int> getDataIDs() const
   {
      return _dataIDs;
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
    * @brief Returns the design specification for the optimization problem.
    *        Information needed to measure the convergence.
    *        In case of manifold mapping it also returns the design specification
    *        for the surrogate model which is updated in every iteration.
    */ // TODO: change to call by ref when Eigen is used.
   virtual std::map<int, utils::DynVector> getDesignSpecification(DataMap& cplData);


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

   typedef tarch::la::DynamicVector<double> DataValues;
   typedef tarch::la::DynamicColumnMatrix<double> DataMatrix;
   typedef tarch::la::DynamicMatrix<double> Matrix;

   /// @brief Logging device.
   static tarch::logging::Log _log;

   /// @brief preconditioner for least-squares system if vectorial system is used.
   PtrPreconditioner _preconditioner;

   /// @brief Cosntant relaxation factor used for first iteration.
   double _initialRelaxation;

   /// @brief Maximum number of old data iterations kept.
   int _maxIterationsUsed;

   /// @brief Maximum number of old timesteps (with data values) kept.
   int _timestepsReused;

   /// @brief Data IDs of data to be involved in the IQN algorithm.
   std::vector<int> _dataIDs;

   /// @brief Data IDs of data not involved in IQN coefficient computation.
   std::vector<int> _secondaryDataIDs;

   /// @brief Indicates the first iteration, where constant relaxation is used.
   bool _firstIteration;

   /* @brief Indicates the first time step, where constant relaxation is used
    *        later, we replace the constant relaxation by a qN-update from last time step.
    */
   bool _firstTimeStep;

   /*
    * @brief If in master-slave mode: True if this process has nodes at the coupling interface
    *        If in server mode: Always true.
    */
   bool _hasNodesOnInterface;

   /* @brief If true, the QN-scheme always performs a underrelaxation in the first iteration of
    *        a new time step. Otherwise, the LS system from the previous time step is used in the
    *        first iteration.
    */
   bool _forceInitialRelaxation;

   /// @brief Solver output from last iteration.
   DataValues _oldXTilde;

   /// @brief Current iteration residuals of IQN data. Temporary.
   DataValues _residuals;

   /// @brief Current iteration residuals of secondary data.
   std::map<int,DataValues> _secondaryResiduals;

   /// @brief Stores residual deltas.
   DataMatrix _matrixV;

   /// @brief Stores x tilde deltas, where x tilde are values computed by solvers.
   DataMatrix _matrixW;
   
   /// @brief Stores the current QR decomposition ov _matrixV, can be updated via deletion/insertion of columns
   QRFactorization _qrV;
   

   /** @brief Indices (of columns in W, V matrices) of 1st iterations of timesteps.
    *
    * When old timesteps are reused (_timestepsReused > 0), the indices of the
    * first iteration of each timestep needs to be stored, such that, e.g., all
    * iterations of the last timestep, or one specific iteration that leads to
    * a singular matrix in the QR decomposition can be removed and tracked.
    */
   std::deque<int> _matrixCols;
   
   /** @brief only needed for the parallel master-slave mode. stores the local dimensions,
    *  i.e., the offsets in _invJacobian for all processors
    */
   std::vector<int> _dimOffsets;


   /** @brief: computes number of cols in least squares system, i.e, number of cols in
    *  _matrixV, _matrixW, _qrV, etc..
    *	 This is necessary only for master-slave mode, when some procs do not have
    *	 any nodes on the coupling interface. In this case, the matrices are not
    *  constructed and we have no information about the number of cols. This info
    *  is needed for master-slave communication. Number of its =! _cols in general.
    */
   int getLSSystemCols();
   int getLSSystemRows();

   /**
     * @brief Marks a iteration sequence as converged.
     *
     * called by the iterationsConverged() method in the BaseQNPostProcessing class
     * handles the postprocessing sepcific action after the convergence of one iteration
     */
    virtual void specializedIterationsConverged(DataMap& cplData) = 0;

   /// @brief updates the V, W matrices (as well as the matrices for the secondary data)
   virtual void updateDifferenceMatrices(DataMap & cplData);
   
   /// @brief concatenates all coupling data involved in the QN system in a single vector
   virtual void concatenateCouplingData(DataMap & cplData);

   /// @brief splits up QN system vector back into the coupling data
   virtual void splitCouplingData(DataMap & cplData);

   /// @brief applies the filter method for the least-squares system, defined in the configuration
   virtual void applyFilter();
   
   /// @brief computes underrelaxation for the secondary data
   virtual void computeUnderrelaxationSecondaryData(DataMap& cplData) = 0;

   /// @brief computes the quasi-Newton update using the specified pp scheme (MVQN, IQNILS)
   virtual void computeQNUpdate(DataMap& cplData, DataValues& xUpdate) = 0;
   
   /// @brief Removes one iteration from V,W matrices and adapts _matrixCols.
   virtual void removeMatrixColumn(int columnIndex);
   
   /// @brief writes info to the _infostream (also in parallel)
   void writeInfo(std::string s, bool allProcs = false);

private:

  /// @brief Concatenation of all coupling data involved in the QN system.
  DataValues _values;

  /// @brief Concatenation of all (old) coupling data involved in the QN system.
  DataValues _oldValues;

  /// @brief Difference between solver input and output from last timestep
  DataValues _oldResiduals;

  /** @brief filter method that is used to maintain good conditioning of the least-squares system
   *        Either of two types: QR1FILTER or QR2Filter
   */
  int _filter;

  /** @brief Determines sensitivity when two matrix columns are considered equal.
   *
   * When during the QR decomposition of the V matrix a pivot element smaller
   * than the singularity limit is found, the matrix is considered to be singular
   * and the corresponding (older) iteration is removed.
   */
  double _singularityLimit;

  /**
    * @brief sets the design specification we want to meet for the objective function,
    *     i. e., we want to solve for argmin_x ||R(x) - q||, with R(x) = H(x) - x
    *     Usually we want to solve for a fixed-point of H, thus solving for argmin_x ||R(x)||
    *     with q=0.
    */
  Eigen::VectorXd _designSpecification;

  /** @brief backup of the V,W and matrixCols data structures. Needed for the skipping of
   *  initial relaxation, if previous time step converged within one iteration i.e., V and W
   *  are empty -- in this case restore V and W with time step t-2.
   */
  DataMatrix _matrixVBackup;
  DataMatrix _matrixWBackup;
  std::deque<int> _matrixColsBackup;

  /// @ brief additional debugging info, is not important for computation:
  int its,tSteps;
  int deletedColumns;

  /// @brief write some debug/post processing info to file
  std::fstream _infostream;


};

}}} // namespace precice, cplscheme, impl

#endif /* PRECICE_CPLSCHEME_BASEQNPOSTPROCESSING_HPP_ */
