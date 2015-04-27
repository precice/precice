// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "DimensionsTest.hpp"
#include "utils/Parallel.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Globals.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::utils::tests::DimensionsTest)

namespace precice {
namespace utils {
namespace tests {

tarch::logging::Log DimensionsTest:: _log ( "precice::utils::tests::DimensionsTest" );

DimensionsTest:: DimensionsTest()
:
  TestCase ( "utils::tests::DimensionsTest" )
{}

void DimensionsTest:: run()
{
  PRECICE_MASTER_ONLY {
    testMethod(testLinearizeDelinearize);
  }
}


void DimensionsTest:: testLinearizeDelinearize()
{
  preciceTrace ( "testLinearizeDelinearize()" );
  { // 2D
    using utils::Vector2D;
    validateEquals (linearize(Vector2D(0.0, 0.0)), 0);
    validateEquals (linearize(Vector2D(1.0, 0.0)), 1);
    validateEquals (linearize(Vector2D(0.0, 1.0)), 2);
    validateEquals (linearize(Vector2D(1.0, 1.0)), 3);

    Vector2D delin0 = delinearize(0, 2);
    validateEquals (delin0(0), -1.0);
    validateEquals (delin0(1), -1.0);

    Vector2D delin1 = delinearize(1, 2);
    validateEquals (delin1(0), 1.0);
    validateEquals (delin1(1), -1.0);

    Vector2D delin2 = delinearize(2, 2);
    validateEquals (delin2(0), -1.0);
    validateEquals (delin2(1), 1.0);

    Vector2D delin3 = delinearize(3, 2);
    validateEquals (delin3(0), 1.0);
    validateEquals (delin3(1), 1.0);
  }

  { // 3D
    using utils::Vector3D;
    validateEquals (linearize(Vector3D(0.0, 0.0, 0.0)), 0);
    validateEquals (linearize(Vector3D(1.0, 0.0, 0.0)), 1);
    validateEquals (linearize(Vector3D(0.0, 1.0, 0.0)), 2);
    validateEquals (linearize(Vector3D(1.0, 1.0, 0.0)), 3);
    validateEquals (linearize(Vector3D(0.0, 0.0, 1.0)), 4);
    validateEquals (linearize(Vector3D(1.0, 0.0, 1.0)), 5);
    validateEquals (linearize(Vector3D(0.0, 1.0, 1.0)), 6);
    validateEquals (linearize(Vector3D(1.0, 1.0, 1.0)), 7);

    Vector3D delin0 = delinearize(0, 3);
    validateEquals (delin0(0), -1.0);
    validateEquals (delin0(1), -1.0);
    validateEquals (delin0(2), -1.0);

    Vector3D delin1 = delinearize(1, 3);
    validateEquals (delin1(0), 1.0);
    validateEquals (delin1(1), -1.0);
    validateEquals (delin1(2), -1.0);

    Vector3D delin2 = delinearize(2, 3);
    validateEquals (delin2(0), -1.0);
    validateEquals (delin2(1), 1.0);
    validateEquals (delin2(2), -1.0);

    Vector3D delin3 = delinearize(3, 3);
    validateEquals (delin3(0), 1.0);
    validateEquals (delin3(1), 1.0);
    validateEquals (delin3(2), -1.0);

    Vector3D delin4 = delinearize(4, 3);
    validateEquals (delin4(0), -1.0);
    validateEquals (delin4(1), -1.0);
    validateEquals (delin4(2), 1.0);

    Vector3D delin5 = delinearize(5, 3);
    validateEquals (delin5(0), 1.0);
    validateEquals (delin5(1), -1.0);
    validateEquals (delin5(2), 1.0);

    Vector3D delin6 = delinearize(6, 3);
    validateEquals (delin6(0), -1.0);
    validateEquals (delin6(1), 1.0);
    validateEquals (delin6(2), 1.0);

    Vector3D delin7 = delinearize(7, 3);
    validateEquals (delin7(0), 1.0);
    validateEquals (delin7(1), 1.0);
    validateEquals (delin7(2), 1.0);
  }
}

//void DimensionsTest:: testGetHyperfaceCornerIndices ()
//{
//  { // 2D
//    int cornerIndices[CompilePower<2,1>::VALUE];
//    getHyperfaceCornerIndices (0, cornerIndices);
//    validateEquals (cornerIndices[0], 0);
//    validateEquals (cornerIndices[1], 2);
//
//    getHyperfaceCornerIndices (1, cornerIndices);
//    validateEquals (cornerIndices[0], 1);
//    validateEquals (cornerIndices[1], 3);
//
//    getHyperfaceCornerIndices (2, cornerIndices);
//    validateEquals (cornerIndices[0], 0);
//    validateEquals (cornerIndices[1], 1);
//
//    getHyperfaceCornerIndices (3, cornerIndices);
//    validateEquals (cornerIndices[0], 2);
//    validateEquals (cornerIndices[1], 3);
//  }
//
//  { // 3D
//    int cornerIndices[CompilePower<2,2>::VALUE];
//    getHyperfaceCornerIndices (0, cornerIndices);
//    validateEquals (cornerIndices[0], 0);
//    validateEquals (cornerIndices[1], 2);
//    validateEquals (cornerIndices[2], 4);
//    validateEquals (cornerIndices[3], 6);
//
//    getHyperfaceCornerIndices (1, cornerIndices);
//    validateEquals (cornerIndices[0], 1);
//    validateEquals (cornerIndices[1], 3);
//    validateEquals (cornerIndices[2], 5);
//    validateEquals (cornerIndices[3], 7);
//
//    getHyperfaceCornerIndices (3, cornerIndices);
//    validateEquals (cornerIndices[0], 2);
//    validateEquals (cornerIndices[1], 3);
//    validateEquals (cornerIndices[2], 6);
//    validateEquals (cornerIndices[3], 7);
//
//    getHyperfaceCornerIndices (5, cornerIndices);
//    validateEquals (cornerIndices[0], 4);
//    validateEquals (cornerIndices[1], 5);
//    validateEquals (cornerIndices[2], 6);
//    validateEquals (cornerIndices[3], 7);
//  }
//}

}}} // namespace precice, utils, tests
