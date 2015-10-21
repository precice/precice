// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_MINITERATIONCONVERGENCEMEASURE_HPP_
#define PRECICE_CPLSCHEME_MINITERATIONCONVERGENCEMEASURE_HPP_

#include "ConvergenceMeasure.hpp"
#include "cplscheme/CouplingData.hpp"
#include "tarch/logging/Log.h"

namespace precice {
namespace cplscheme {
namespace impl {

class MinIterationConvergenceMeasure : public ConvergenceMeasure
{
public:

   MinIterationConvergenceMeasure ( int minimumIterationCount );

   virtual ~MinIterationConvergenceMeasure() {}

   virtual void newMeasurementSeries();

   virtual void measure (
      const utils::DynVector & oldValues,
      const utils::DynVector & newValues,
      const utils::DynVector& designSpecification)
   {
     preciceTrace("measure()");
     _currentIteration++;
     _isConvergence = _minimumIterationCount <= _currentIteration
                      ? true
                      : false;
     preciceDebug("Iteration number = " << _currentIteration
                  << ", convergence = " << _isConvergence);
   }

   virtual bool isConvergence() const
   {
     return _isConvergence;
   }

   virtual std::string printState()
   {
     std::ostringstream os;
     os << "min iteration convergence measure: ";
     os << "#it = " << _currentIteration << " of " << _minimumIterationCount;
     os << ", conv = ";
     if (_isConvergence) os << "true";
     else os << "false";
     return os.str();
   }

private:

  static tarch::logging::Log _log;

  int _minimumIterationCount;

  int _currentIteration;

  bool _isConvergence;
};

}}} // namespace precice, cplscheme, impl

#endif // PRECICE_CPLSCHEME_MINITERATIONCONVERGENCEMEASURE_HPP_
