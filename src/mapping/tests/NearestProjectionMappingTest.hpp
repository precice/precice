// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_MAPPING_NEARESTPROJETIONMAPPINGTEST_HPP_
#define PRECICE_MAPPING_NEARESTPROJETIONMAPPINGTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace mapping {
namespace tests {


class NearestProjectionMappingTest
:
   public tarch::tests::TestCase
{
public:

     NearestProjectionMappingTest ();

   virtual ~NearestProjectionMappingTest() {};

   virtual void setUp () {}

   virtual void run ();

private:

   static tarch::logging::Log _log;

   void testConservativeNonIncremental();

   void testConsistentNonIncremental();
};

}}} // namespace precice, mapping, tests

#endif /* PRECICE_MAPPING_NEARESTPROJETIONMAPPINGTEST_HPP_ */
