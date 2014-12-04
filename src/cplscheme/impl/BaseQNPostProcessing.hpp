// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_BASEQNPOSTPROCESSING_HPP_
#define PRECICE_CPLSCHEME_BASEQNPOSTPROCESSING_HPP_

#include "PostProcessing.hpp"
#include "mesh/SharedPointer.hpp"
#include "tarch/logging/Log.h"
#include "tarch/la/DynamicColumnMatrix.h"
#include "tarch/la/DynamicVector.h"
#include <deque>

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
      int    maxIterationsUsed,
      int    timestepsReused,
      double singularityLimit,
      std::vector<int>    dataIDs,
      std::map<int,double>    scalings);

   /**
    * @brief Destructor, empty.
    */
   virtual ~BaseQNPostProcessing() {}

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
   
   virtual void performPPSecondaryData(DataMap& cplData) = 0;

   /**
    * @brief Marks a iteration sequence as converged.
    *
    * Since convergence measurements are done outside the post-processing, this
    * method has to be used to signalize convergence to the post-processing.
    */
   virtual void iterationsConverged(DataMap& cplData);

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
   

protected:

   typedef tarch::la::DynamicVector<double> DataValues;

   typedef tarch::la::DynamicColumnMatrix<double> DataMatrix;

   // @brief Logging device.
   static tarch::logging::Log _log;

   // @brief Cosntant relaxation factor used for first iteration.
   double _initialRelaxation;

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

   // @brief Data IDs of data to be involved in the IQN algorithm.
   std::vector<int> _dataIDs;

   // @brief Data IDs of data not involved in IQN coefficient computation.
   std::vector<int> _secondaryDataIDs;

   // @brief Scales data by fixed value.
   //
   // When more than one data is used to compute the IQN linear combination of
   // old data columns (_dataIDs.size() > 0), all data should have similar
   // magnitude, in order to be similarly important in the least-squares solution.
   std::map<int,double> _scalings;

   // @brief Indicates the first iteration, where constant relaxation is used.
   bool _firstIteration;
   bool _firstTimeStep;

   // @brief Solver output from last iteration.
   DataValues _oldXTilde;

   // @brief Secondary data solver output from last iteration.
   //std::map<int,DataValues> _secondaryOldXTildes;

   // @brief Current iteration residuals of IQN data. Temporary.
   DataValues _residuals;

   // @brief Current iteration residuals of secondary data.
   std::map<int,DataValues> _secondaryResiduals;

   // @brief Temporary used in performPostProcessing().
   DataValues _scaledValues;

   // @brief Temporary used in performPostProcessing().
   DataValues _scaledOldValues;

   // @brief Difference between solver input and output from last timestep
   DataValues _oldResiduals;

   // @brief Stores residual deltas.
   DataMatrix _matrixV;

   // @brief Stores x tilde deltas, where x tilde are values computed by solvers.
   DataMatrix _matrixW;
   
   DataMatrix _matrixVBackup;
   DataMatrix _matrixWBackup;
   std::deque<int> _matrixColsBackup;

   // @brief Indices (of columns in W, V matrices) of 1st iterations of timesteps.
   //
   // When old timesteps are reused (_timestepsReused > 0), the indices of the
   // first iteration of each timestep needs to be stored, such that, e.g., all
   // iterations of the last timestep, or one specific iteration that leads to
   // a singular matrix in the QR decomposition can be removed and tracked.
   std::deque<int> _matrixCols;

   // @brief computes the quasi-Newton update using the specified pp scheme (MVQN, IQNILS)
   virtual void computeQNUpdate(DataMap& cplData, DataValues& xUpdate) = 0;
   
   // @brief Removes one iteration from V,W matrices and adapts _matrixCols.
   virtual void removeMatrixColumn(int columnIndex);
   
};

}}} // namespace precice, cplscheme, impl

#endif /* PRECICE_CPLSCHEME_BASEQNPOSTPROCESSING_HPP_ */
