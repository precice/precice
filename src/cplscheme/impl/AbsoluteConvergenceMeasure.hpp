// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef ABSOLUTECONVERGENCEMEASURE_HPP_
#define ABSOLUTECONVERGENCEMEASURE_HPP_

#include "ConvergenceMeasure.hpp"
#include "utils/Helpers.hpp"
#include "utils/Dimensions.hpp"
//#include "utils/NumericalCompare.hpp"
#include "tarch/logging/Log.h"
#include "utils/MasterSlave.hpp"

namespace precice {
   namespace cplscheme {
      namespace tests {
         class AbsoluteConvergenceMeasureTest;
      }
   }
}

// ------------------------------------------------------------ CLASS DEFINTION

namespace precice {
namespace cplscheme {
namespace impl {

/**
 * @brief Measures the convergence from an old data set to a new one.
 *
 * The convergence is evaluated by looking at the two norm of the differences
 * between each data value from the new and old data set. If the two norm is
 * equal or below a given limit, convergence is achieved.
 *
 * For a description of how to perform the measurement, see class
 * ConvergenceMeasure.
 */
class AbsoluteConvergenceMeasure : public ConvergenceMeasure
{
public:

   AbsoluteConvergenceMeasure ( double convergenceLimit );

   virtual ~AbsoluteConvergenceMeasure() {};

   virtual void newMeasurementSeries ()
   {
      _isConvergence = false;
   }

   virtual void measure (
      const utils::DynVector& oldValues,
      const utils::DynVector& newValues )
   {
      _normDiff = utils::MasterSlave::l2norm(newValues - oldValues);
      _isConvergence = _normDiff <= _convergenceLimit;
//      preciceInfo ( "measure()", "Absolute convergence measure: "
//                     << "two-norm differences = " << normDiff
//                     << ", convergence limit = " << _convergenceLimit
//                     << ", convergence = " << _isConvergence );
   }

   virtual bool isConvergence () const
   {
      return _isConvergence;
   }

   /**
    * @brief Adds current convergence information to output stream.
    */
   virtual std::string printState()
   {
     std::ostringstream os;
     os << "absolute convergence measure: ";
     os << "two-norm diff = " << _normDiff;
     os << ", limit = " << _convergenceLimit;
     os << ", conv = ";
     if (_isConvergence) os << "true";
     else os << "false";
     return os.str();
   }
   
   virtual double getNormResidual()
   {
    return _normDiff; 
   }

private:

   static tarch::logging::Log _log;

   double _convergenceLimit;

   double _normDiff;

   bool _isConvergence;
};

}}} // namespace precice, cplscheme, impl

#endif /* ABSOLUTECONVERGENCEMEASURE_HPP_ */
