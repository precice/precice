// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_TESTS_RELATIVECONVERGENCEMEASURETEST_HPP_
#define PRECICE_CPLSCHEME_TESTS_RELATIVECONVERGENCEMEASURETEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace cplscheme {
namespace tests {

/**
 * @brief Tests class RelativeConvergenceMeasure.
 */
class RelativeConvergenceMeasureTest : public tarch::tests::TestCase
{
public:

   /**
    * @brief Constructor.
    */
   RelativeConvergenceMeasureTest ();

   /**
    * @brief Destructor.
    */
   virtual ~RelativeConvergenceMeasureTest () {};

   /**
    * @brief Empty.
    */
   virtual void setUp () {}

   /**
    * @brief Calls all test methods.
    */
   virtual void run ();

private:

   // @brief Logging device.
   static tarch::logging::Log _log;

   /**
    * @brief Tests convergence measurement with Vector datasets.
    */
   void testMeasureData ();

   /**
    * @brief Tests convergence measurement with double datasets.
    */
//   void testMeasureDoubleData ();

   /**
    * @brief Tests convergence measurement with int datasets.
    */
//   void testMeasureIntegerData ();
};

}}} // namespace precice, cplscheme, tests

#endif // PRECICE_CPLSCHEME_TESTS_RELATIVECONVERGENCEMEASURETEST_HPP_
