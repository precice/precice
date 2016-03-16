// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_GEOMETRY_TESTS_IMPORTGEOMETRYTEST_HPP_
#define PRECICE_GEOMETRY_TESTS_IMPORTGEOMETRYTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace geometry {
namespace tests {

/**
 * @brief Provides tests for class ImportGeometry.
 */
class ImportGeometryTest : public tarch::tests::TestCase
{
public:

   ImportGeometryTest ();

   virtual ~ImportGeometryTest() {};

   /**
    * @brief Empty.
    */
   virtual void setUp () {}

   virtual void run ();

private:

   static logging::Logger _log;

   void testImportVRMLConfig ();
};

}}} // namespace precice, geometry, tests

#endif /* PRECICE_GEOMETRY_TESTS_IMPORTGEOMETRYTEST_HPP_ */
