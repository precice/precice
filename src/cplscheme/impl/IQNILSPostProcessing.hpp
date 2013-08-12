// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_IQNILSPOSTPROCESSING_HPP_
#define PRECICE_CPLSCHEME_IQNILSPOSTPROCESSING_HPP_

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

class IQNILSPostProcessing : public PostProcessing
{
public:

   IQNILSPostProcessing (
      double initialRelaxation,
      int    maxIterationsUsed,
      int    timestepsReused,
      double singularityLimit,
      std::vector<int>    dataIDs );

   virtual ~IQNILSPostProcessing() {};

   virtual std::vector<int> getDataIDs () const
   {
      return _dataIDs;
   }

   virtual void initialize (
      DataMap & cplData );

   virtual void performPostProcessing (
      DataMap & cplData );

   virtual void iterationsConverged (
      DataMap & cplData );

   virtual void exportState(io::TXTWriter& writer);

   virtual void importState(io::TXTReader& reader);

private:

   typedef tarch::la::DynamicVector<double> DataValues;

   typedef tarch::la::DynamicColumnMatrix<double> DataMatrix;

   static tarch::logging::Log _log;

   double _initialRelaxation;

   int _maxIterationsUsed;

   int _timestepsReused;

   double _singularityLimit;

   std::vector<int> _dataIDs;

   bool _firstIteration;

   // @brief Solver output from last iteration.
   DataValues _oldXTilde;

   // @brief Difference between solver input and output from last timestep
   DataValues _oldResiduals;

   // @brief Stores residual deltas.
   DataMatrix _matrixV;

   // @brief Stores x tilde deltas, where x tilde are values computed by solvers.
   DataMatrix _matrixW;

   std::deque<int> _matrixCols;

   void removeMatrixColumn ( int columnIndex );
};

}}} // namespace precice, cplscheme, impl

#endif /* PRECICE_CPLSCHEME_IQNILSPOSTPROCESSING_HPP_ */
