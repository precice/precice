// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "PointerVectorTest.hpp"
#include "../PointerVector.hpp"
#include "../Parallel.hpp"
#include "../Helpers.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::utils::tests::PointerVectorTest)

namespace precice {
namespace utils {
namespace tests {

tarch::logging::Log PointerVectorTest:: _log ( "precice::utils::tests::PointerVectoTest" );

PointerVectorTest:: PointerVectorTest ()
:
   TestCase ( "utils::tests::PointerVectorTest" )
{}

void PointerVectorTest:: run ()
{
   PRECICE_MASTER_ONLY {
      preciceTrace ( "run()" );
      ptr_vector<double> ptrVector;

   }
}

}}} // namespace precice, utils, tests
