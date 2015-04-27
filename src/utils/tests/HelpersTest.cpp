// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "HelpersTest.hpp"
#include "utils/Helpers.hpp"
#include "utils/Parallel.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Globals.hpp"
#include <vector>
#include <list>
#include <map>

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::utils::tests::HelpersTest)

namespace precice {
namespace utils {
namespace tests {

tarch::logging::Log HelpersTest:: _log ("precice::utils::HelpersTest");

HelpersTest:: HelpersTest()
:
  TestCase ("utils::HelpersTest")
{}

void HelpersTest:: run()
{
  PRECICE_MASTER_ONLY {
    testMethod(testAppendTo);
    testMethod(testOperatorPlusForVectors);
  }
}

void HelpersTest:: testAppendTo()
{
  preciceTrace ( "testAppendTo()" );
  // Test std::vector
  std::vector<double> doubleVector;
  appendTo(doubleVector) ( 1.0 )( 2.0 )( 3.0 );
  validateEquals ( doubleVector.size(), 3 );
  validateNumericalEquals ( doubleVector[0], 1.0 );
  validateNumericalEquals ( doubleVector[1], 2.0 );
  validateNumericalEquals ( doubleVector[2], 3.0 );
  appendTo(doubleVector) ( 4.0 );
  validateEquals ( doubleVector.size(), 4 );
  validateNumericalEquals ( doubleVector[3], 4.0 );

  std::vector<int> intVector;
  appendTo(intVector) ( 1 )( 2 )( 3 );
  validateEquals ( intVector.size(), 3 );
  validateNumericalEquals ( intVector[0], 1 );
  validateNumericalEquals ( intVector[1], 2 );
  validateNumericalEquals ( intVector[2], 3 );
  appendTo(intVector) ( 4 );
  validateEquals ( intVector.size(), 4 );
  validateNumericalEquals ( intVector[3], 4 );

  using utils::Vector2D;
  std::vector<Vector2D> vectorVector;
  appendTo(vectorVector) ( Vector2D(1.0) )( Vector2D(2.0) )( Vector2D(3.0) );
  validateEquals ( vectorVector.size(), 3 );
  validate ( tarch::la::equals(vectorVector[0], Vector2D(1.0)) );
  validate ( tarch::la::equals(vectorVector[1], Vector2D(2.0)) );
  validate ( tarch::la::equals(vectorVector[2], Vector2D(3.0)) );
  appendTo(vectorVector) ( Vector2D(4.0) );
  validateEquals ( vectorVector.size(), 4 );
  validate ( tarch::la::equals(vectorVector[3], Vector2D(4.0)) );

  // Test std::list
  std::list<double> doubleList;
  appendTo(doubleList) ( 1.0 )( 2.0 )( 3.0 );
  validateEquals ( doubleList.size(), 3 );
  appendTo(doubleList) ( 4.0 );
  validateEquals ( doubleList.size(), 4 );
  std::list<double>::iterator doubleListIter = doubleList.begin ();
  validateNumericalEquals ( *doubleListIter, 1.0 );
  doubleListIter++;
  validateNumericalEquals ( *doubleListIter, 2.0 );
  doubleListIter++;
  validateNumericalEquals ( *doubleListIter, 3.0 );
  doubleListIter++;
  validateNumericalEquals ( *doubleListIter, 4.0 );

  std::list<int> intList;
  appendTo(intList) ( 1 )( 2 )( 3 );
  validateEquals ( intList.size(), 3 );
  appendTo(intList) ( 4 );
  validateEquals ( intList.size(), 4 );
  std::list<int>::iterator intListIter = intList.begin ();
  validateNumericalEquals ( *intListIter, 1 );
  intListIter++;
  validateNumericalEquals ( *intListIter, 2 );
  intListIter++;
  validateNumericalEquals ( *intListIter, 3 );
  intListIter++;
  validateNumericalEquals ( *intListIter, 4 );

  std::list<Vector2D> vectorList;
  appendTo(vectorList) ( Vector2D(1.0) )( Vector2D(2.0) )( Vector2D(3.0) );
  validateEquals ( vectorList.size(), 3 );
  appendTo(vectorList) ( Vector2D(4.0) );
  validateEquals ( vectorList.size(), 4 );
  std::list<Vector2D>::iterator vectorListIter = vectorList.begin ();
  validate ( tarch::la::equals(*vectorListIter, Vector2D(1.0)) );
  vectorListIter++;
  validate ( tarch::la::equals(*vectorListIter, Vector2D(2.0)) );
  vectorListIter++;
  validate ( tarch::la::equals(*vectorListIter, Vector2D(3.0)) );
  vectorListIter++;
  validate ( tarch::la::equals(*vectorListIter, Vector2D(4.0)) );

  // Test std::map
  std::map<int,Vector2D> aMap;
  appendTo(aMap) ( std::make_pair<int,Vector2D>(1, Vector2D(1.0)) )
                    ( std::make_pair<int,Vector2D>(2, Vector2D(2.0)) )
                    ( std::make_pair<int,Vector2D>(3, Vector2D(3.0))  );
  validateEquals ( aMap.size(), 3 );
  validate ( aMap.find(1) != aMap.end() );
  validate ( aMap.find(2) != aMap.end() );
  validate ( aMap.find(3) != aMap.end() );
  validate ( tarch::la::equals(aMap[1], Vector2D(1.0)) );
  validate ( tarch::la::equals(aMap[2], Vector2D(2.0)) );
  validate ( tarch::la::equals(aMap[3], Vector2D(3.0)) );
}

void HelpersTest:: testOperatorPlusForVectors()
{
  preciceTrace ( "testOperatorPlusForVectors()" );
  std::vector<double> doubleVector;
  doubleVector += 1.0, 2.0, 3.0;
  validateEquals ( doubleVector.size(), 3 );
  validateEquals ( doubleVector[0], 1.0 );
  validateEquals ( doubleVector[1], 2.0 );
  validateEquals ( doubleVector[2], 3.0 );

  using utils::Vector2D;
  std::vector<Vector2D> vectorVector;
  vectorVector += Vector2D(1.0), Vector2D(2.0), Vector2D(3.0);
  validateEquals ( vectorVector.size(), 3 );
  validate ( tarch::la::equals(vectorVector[0], Vector2D(1.0)) );
  validate ( tarch::la::equals(vectorVector[1], Vector2D(2.0)) );
  validate ( tarch::la::equals(vectorVector[2], Vector2D(3.0)) );
}

}}} // namespace precice, utils, tests
