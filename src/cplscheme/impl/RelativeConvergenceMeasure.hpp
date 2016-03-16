// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef RELATIVECONVERGENCEMEASURE_HPP_
#define RELATIVECONVERGENCEMEASURE_HPP_

#include "ConvergenceMeasure.hpp"
#include "../CouplingData.hpp"
#include "utils/Helpers.hpp"
#include "utils/Dimensions.hpp"
//#include "utils/NumericalCompare.hpp"
#include "logging/Logger.hpp"
#include "utils/MasterSlave.hpp"
#include "tarch/la/ScalarOperations.h"

namespace precice {
   namespace cplscheme {
      namespace tests {
         class RelativeConvergenceMeasureTest;
      }
   }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace cplscheme {
namespace impl {

/**
 * @brief Measures the convergence from an old data set to a new one.
 *
 * The convergence is evaluated by looking at the two norm of the differences
 * between each data value from the new and old data set. If the two norm is
 * equal or below a given percentage of the norm of the old data set,
 * convergence is achieved.
 *
 * For a description of how to perform the measurement, see class
 * ConvergenceMeasure.
 */
class RelativeConvergenceMeasure : public ConvergenceMeasure
{
public:

   /**
    * @brief Constructor.
    *
    * @param convergenceLimitPercent [IN]
    *        Limit to define convergence relative to the norm of the current
    *        new dataset. Has to be in $] 0 ; 1 ]$.
    */
   RelativeConvergenceMeasure ( double convergenceLimitPercent );

   virtual ~RelativeConvergenceMeasure () {};

   virtual void newMeasurementSeries ()
   {
      _isConvergence = false;
   }

   virtual void measure (
      const  Eigen::VectorXd& oldValues,
      const  Eigen::VectorXd& newValues,
      const  Eigen::VectorXd& designSpecification)
   {
/*
     std::cout<<"\n-------"<<std::endl;
     std::cout<<"   old val: \n"<<oldValues<<std::endl;
     std::cout<<"   new val: \n"<<newValues<<"\n"<<std::endl;
     std::cout<<"   design spec: \n"<<designSpecification<<"\n"<<std::endl;
     std::cout<<"-------\n"<<std::endl;
*/

     _normDiff = utils::MasterSlave::l2norm((newValues - oldValues) - designSpecification);
     _norm = utils::MasterSlave::l2norm(newValues + designSpecification);
     _isConvergence = _normDiff <= _norm * _convergenceLimitPercent;
//      preciceInfo ( "measure()", "Relative convergence measure: "
//                    << "two-norm differences = " << normDiff
//                    << ", convergence limit = "
//                    << normNew * _convergenceLimitPercent
//                    << ", convergence = " << _isConvergence );
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
     os << "relative convergence measure: ";
     os << "two-norm diff = " << _normDiff;
     os << ", limit = " << _norm * _convergenceLimitPercent;
     os << ", conv = ";
     if (_isConvergence) os << "true";
     else os << "false";
     return os.str();
   }
   
   virtual double getNormResidual()
   {
	   if (tarch::la::equals(_norm, 0.))
		   return std::numeric_limits<double>::infinity();
	   else
		   return _normDiff/_norm;
   }

private:

   static logging::Logger _log;

   double _convergenceLimitPercent;

   double _normDiff;

   double _norm;

   bool _isConvergence;
};

}}} // namespace precice, cplscheme, impl

#endif /* RELATIVECONVERGENCEMEASURE_HPP_ */
