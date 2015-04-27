// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_ITESTS_SCALEACTIONTEST_HPP_
#define PRECICE_ITESTS_SCALEACTIONTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace action {
namespace tests {


/**
 * @brief Unit test class for class action::ScaleProperty.
 */
class ScaleActionTest : public tarch::tests::TestCase
{
public:

   /**
    * @brief Constructor.
    */
   ScaleActionTest ();

   /**
    * @brief Destructor.
    */
   virtual ~ScaleActionTest () {};

   /**
    * @brief Empty.
    */
   virtual void setUp () {}

   /**
    * @brief All test methods are called from here.
    */
   virtual void run ();

private:

   // @brief Logging device.
   static tarch::logging::Log _log;

   /**
    * @brief Test division of data by area of mesh elements (edges or triangles).
    */
   void testDivideByArea ();

   void testScaleByComputedTimestepLength ();

   /**
    * @brief Tests scaling by ratio of computed dt to full dt length.
    */
   void testScaleByComputedTimestepPartLength ();

   /**
    * @brief Tests configuration from xml files.
    */
   void testConfiguration ();
};

}}} // namespace precice, action, tests

#endif /* PRECICE_ITESTS_SCALEACTIONTEST_HPP_ */
