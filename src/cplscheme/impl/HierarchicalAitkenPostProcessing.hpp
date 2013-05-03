// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_HIERARCHICALAITKENPOSTPROCESSING_HPP_
#define PRECICE_CPLSCHEME_HIERARCHICALAITKENPOSTPROCESSING_HPP_

#include "PostProcessing.hpp"
#include "tarch/logging/Log.h"
#include "tarch/la/DynamicVector.h"
#include <vector>

namespace precice {
namespace cplscheme {
namespace impl {

class HierarchicalAitkenPostProcessing : public PostProcessing
{
public:

  HierarchicalAitkenPostProcessing (
    double initialRelaxation,
    int    dataID );

  virtual ~HierarchicalAitkenPostProcessing () {};

  virtual int getDataID () const
  {
    return _dataID;
  }

  virtual void initialize ( DataMap & cplData );

  virtual void performPostProcessing ( DataMap & cplData );

  virtual void iterationsConverged ( DataMap & cplData);

private:

  static tarch::logging::Log _log;

  double _initialRelaxation;

  int _dataID;

  std::vector<double> _aitkenFactors;

  int _iterationCounter;

  tarch::la::DynamicVector<double> _residual;

  void computeAitkenFactor (
    size_t level,
    double nominator,
    double denominator );
};

}}} // namespace precice, cplscheme, impl

#endif /* PRECICE_CPLSCHEME_HIERARCHICALAITKENPOSTPROCESSING_HPP_ */
