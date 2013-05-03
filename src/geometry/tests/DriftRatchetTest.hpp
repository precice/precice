// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_GEOMETRY_TESTS_DRIFTRATCHETTEST_HPP_
#define PRECICE_GEOMETRY_TESTS_DRIFTRATCHETTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace geometry {
namespace tests {


class DriftRatchetTest : public tarch::tests::TestCase
{
public:

   static tarch::logging::Log _log;

   DriftRatchetTest ();

   virtual ~DriftRatchetTest () {};

   /**
    * @brief Empty.
    */
   virtual void setUp () {}

   virtual void run ();
};

}}} // namespace precice, geometry, tests

#endif // PRECICE_GEOMETRY_TESTS_DRIFTRATCHETTEST_HPP_
