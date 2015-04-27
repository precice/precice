// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "CompositionalCouplingSchemeTest.hpp"
#include "DummyCouplingScheme.hpp"
#include "../SharedPointer.hpp"
#include "../CompositionalCouplingScheme.hpp"
#include "../Constants.hpp"
#include "../config/CouplingSchemeConfiguration.hpp"
#include "mesh/PropertyContainer.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "geometry/config/GeometryConfiguration.hpp"
#include "com/Communication.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "m2n/M2N.hpp"
#include "utils/xml/XMLTag.hpp"
#include "utils/Dimensions.hpp"
#include <vector>

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::cplscheme::tests::CompositionalCouplingSchemeTest)

namespace precice {
namespace cplscheme {
namespace tests {

using utils::Vector3D;

tarch::logging::Log CompositionalCouplingSchemeTest::
   _log ( "precice::cplscheme::tests::CompositionalCouplingSchemeTest" );

CompositionalCouplingSchemeTest:: CompositionalCouplingSchemeTest()
:
   TestCase("cplscheme::CompositionalCouplingSchemeTest"),
   _pathToTests()
{}

void CompositionalCouplingSchemeTest:: setUp ()
{
  _pathToTests = utils::Globals::getPathToSources() + "/cplscheme/tests/";
}

void CompositionalCouplingSchemeTest:: run ()
{
# ifndef PRECICE_NO_MPI
  typedef utils::Parallel Par;
  PRECICE_MASTER_ONLY{
    testMethod(testDummySchemeCompositions);
  }
  if (Par::getCommunicatorSize() > 2){
    std::vector<int> ranks;
    ranks += 0, 1, 2;
    Par::Communicator comm = Par::getRestrictedCommunicator(ranks);
    if (Par::getProcessRank() <= 2){
      Par::setGlobalCommunicator(comm) ;
      validateEquals(Par::getCommunicatorSize(), 3);
      testMethod(testExplicitSchemeComposition1);
      //TODO did produce a deadlock on Benjamin's laptop
      //testMethod(testImplicitSchemeComposition);
      //testMethod(testImplicitExplicitSchemeComposition);
      //testMethod(testExplicitImplicitSchemeComposition);
      Par::setGlobalCommunicator(Par::getCommunicatorWorld());
    }
  }
# endif // not PRECICE_NO_MPI
}

void CompositionalCouplingSchemeTest:: testDummySchemeCompositions()
{
  preciceTrace("testDummySchemeCompositions()");
  std::string writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());
  std::string readIterationCheckpoint(constants::actionReadIterationCheckpoint());
  { // Test one explicit dummy coupling scheme
    preciceDebug("Test E");
    int numberIterations = 1;
    int maxTimesteps = 10;
    PtrCouplingScheme scheme(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    CompositionalCouplingScheme composition;
    composition.addCouplingScheme(scheme);
    composition.initialize(0.0, 1);
    int advances = 0;
    while (composition.isCouplingOngoing()){
      composition.advance();
      advances++;
    }
    composition.finalize();
    validateEquals(advances, 10);
    validateEquals(scheme->getTimesteps()-1, 10);
  }
  { // Test one implicit dummy coupling scheme
    preciceDebug("Test I(2)");
    int numberIterations = 2;
    int maxTimesteps = 10;
    PtrCouplingScheme scheme(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    CompositionalCouplingScheme composition;
    composition.addCouplingScheme(scheme);
    composition.initialize(0.0, 1);
    int advances = 0;
    while (composition.isCouplingOngoing()){
      composition.advance();
      advances++;
    }
    composition.finalize();
    validateEquals(advances, 20);
    validateEquals(scheme->getTimesteps()-1, 10);
  }
  { // Test two explicit dummy coupling schemes
    preciceDebug("Test E, E");
    int numberIterations = 1;
    int maxTimesteps = 10;
    PtrCouplingScheme scheme1(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    PtrCouplingScheme scheme2(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    CompositionalCouplingScheme composition;
    composition.addCouplingScheme(scheme1);
    composition.addCouplingScheme(scheme2);
    composition.initialize(0.0, 1);
    int advances = 0;
    while (composition.isCouplingOngoing()){
      composition.advance();
      advances++;
    }
    composition.finalize();
    validateEquals(advances, 10);
    validateEquals(scheme1->getTimesteps()-1, 10);
    validateEquals(scheme2->getTimesteps()-1, 10);
  }
  { // Test three explicit dummy coupling schemes
    preciceDebug("Test E, E, E");
    int numberIterations = 1;
    int maxTimesteps = 10;
    PtrCouplingScheme scheme1(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    PtrCouplingScheme scheme2(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    PtrCouplingScheme scheme3(
        new DummyCouplingScheme(numberIterations, maxTimesteps));
    CompositionalCouplingScheme composition;
    composition.addCouplingScheme(scheme1);
    composition.addCouplingScheme(scheme2);
    composition.addCouplingScheme(scheme3);
    composition.initialize(0.0, 1);
    int advances = 0;
    while (composition.isCouplingOngoing()){
      composition.advance();
      advances++;
    }
    composition.finalize();
    validateEquals(advances, 10);
    validateEquals(scheme1->getTimesteps()-1, 10);
    validateEquals(scheme2->getTimesteps()-1, 10);
    validateEquals(scheme3->getTimesteps()-1, 10);
  }
  { // Test two implicit dummy coupling schemes
    preciceDebug("Test I(2), I(2)");
    int numberIterations = 2;
    int maxTimesteps = 10;
    PtrCouplingScheme scheme1(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    PtrCouplingScheme scheme2(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    CompositionalCouplingScheme composition;
    composition.addCouplingScheme(scheme1);
    composition.addCouplingScheme(scheme2);
    composition.initialize(0.0, 1);
    int advances = 0;
    while (composition.isCouplingOngoing()){
      composition.advance();
      advances++;
      if (advances%2 == 1){
        validate(scheme1->isActionRequired(readIterationCheckpoint));
        validate(scheme2->isActionRequired(readIterationCheckpoint));
      }
      else if (advances%2 == 0){
        validate(scheme1->isActionRequired(writeIterationCheckpoint));
        validate(scheme2->isActionRequired(writeIterationCheckpoint));
      }
    }
    composition.finalize();
    validateEquals(advances, 20);
    validateEquals(scheme1->getTimesteps()-1, 10);
    validateEquals(scheme2->getTimesteps()-1, 10);
  }
  { // Test two implicit dummy coupling schemes with different iteration number
    preciceDebug("Test I(2), I(3)");
    int numberIterations = 2;
    int maxTimesteps = 10;
    PtrCouplingScheme scheme1(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    numberIterations = 3;
    PtrCouplingScheme scheme2(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    CompositionalCouplingScheme composition;
    composition.addCouplingScheme(scheme1);
    composition.addCouplingScheme(scheme2);
    composition.initialize(0.0, 1);
    int advances = 0;
    while (composition.isCouplingOngoing()){
      composition.advance();
      advances++;
      if (advances%3 == 2){
        validate(scheme1->isActionRequired(writeIterationCheckpoint));
        validate(scheme2->isActionRequired(readIterationCheckpoint));
      }
      else if (advances%3 == 1){
        validate(scheme1->isActionRequired(readIterationCheckpoint));
        validate(scheme2->isActionRequired(readIterationCheckpoint));
      }
      else if (advances%3 == 0){
        validate(scheme1->isActionRequired(writeIterationCheckpoint));
        validate(scheme2->isActionRequired(writeIterationCheckpoint));
      }
    }
    composition.finalize();
    validateEquals(advances, 30);
    validateEquals(scheme1->getTimesteps()-1, 10);
    validateEquals(scheme2->getTimesteps()-1, 10);
  }
  { // Test three implicit dummy coupling schemes
    preciceDebug("Test I(2), I(2), I(2)");
    int numberIterations = 2;
    int maxTimesteps = 10;
    PtrCouplingScheme scheme1(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    PtrCouplingScheme scheme2(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    PtrCouplingScheme scheme3(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    CompositionalCouplingScheme composition;
    composition.addCouplingScheme(scheme1);
    composition.addCouplingScheme(scheme2);
    composition.addCouplingScheme(scheme3);
    composition.initialize(0.0, 1);
    int advances = 0;
    while (composition.isCouplingOngoing()){
      composition.advance();
      advances++;
      if (advances%2 == 0){
        validate(scheme1->isActionRequired(writeIterationCheckpoint));
        validate(scheme2->isActionRequired(writeIterationCheckpoint));
        validate(scheme3->isActionRequired(writeIterationCheckpoint));
      }
      else {
        validate(scheme1->isActionRequired(readIterationCheckpoint));
        validate(scheme2->isActionRequired(readIterationCheckpoint));
        validate(scheme3->isActionRequired(readIterationCheckpoint));
      }
    }
    composition.finalize();
    validateEquals(advances, 20);
    validateEquals(scheme1->getTimesteps()-1, 10);
    validateEquals(scheme2->getTimesteps()-1, 10);
    validateEquals(scheme3->getTimesteps()-1, 10);
  }
  { // Test three implicit dummy coupling schemes
    preciceDebug("Test I(3), I(4), I(2)");
    int numberIterations = 3;
    int maxTimesteps = 10;
    PtrCouplingScheme scheme1(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    numberIterations = 4;
    PtrCouplingScheme scheme2(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    numberIterations = 2;
    PtrCouplingScheme scheme3(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    CompositionalCouplingScheme composition;
    composition.addCouplingScheme(scheme1);
    composition.addCouplingScheme(scheme2);
    composition.addCouplingScheme(scheme3);
    composition.initialize(0.0, 1);
    int advances = 0;
    while (composition.isCouplingOngoing()){
      composition.advance();
      advances++;
      if (advances%4 == 0){
        validate(scheme1->isActionRequired(writeIterationCheckpoint));
        validate(scheme2->isActionRequired(writeIterationCheckpoint));
        validate(scheme3->isActionRequired(writeIterationCheckpoint));
      }
      else if (advances%4 == 1){
        validate(scheme1->isActionRequired(readIterationCheckpoint));
        validate(scheme2->isActionRequired(readIterationCheckpoint));
        validate(scheme3->isActionRequired(readIterationCheckpoint));
      }
      else if (advances%4 == 2){
        validate(scheme1->isActionRequired(readIterationCheckpoint));
        validate(scheme2->isActionRequired(readIterationCheckpoint));
        validate(scheme3->isActionRequired(writeIterationCheckpoint));
      }
      else if (advances%4 == 3){
        validate(scheme1->isActionRequired(writeIterationCheckpoint));
        validate(scheme2->isActionRequired(readIterationCheckpoint));
        validate(scheme3->isActionRequired(writeIterationCheckpoint));
      }
    }
    composition.finalize();
    validateEquals(advances, 40);
    validateEquals(scheme1->getTimesteps()-1, 10);
    validateEquals(scheme2->getTimesteps()-1, 10);
    validateEquals(scheme3->getTimesteps()-1, 10);
  }
  { // Test E, I(2)
    preciceDebug("Test E, I(2)");
    int numberIterations = 1;
    int maxTimesteps = 10;
    PtrCouplingScheme scheme1(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    numberIterations = 2;
    PtrCouplingScheme scheme2(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    CompositionalCouplingScheme composition;
    composition.addCouplingScheme(scheme1);
    composition.addCouplingScheme(scheme2);
    composition.initialize(0.0, 1);
    int advances = 0;
    while (composition.isCouplingOngoing()){
      composition.advance();
      advances++;
      if (advances%2 == 0){
        validate(scheme2->isActionRequired(writeIterationCheckpoint));
        validateEquals(scheme1->getTimesteps()-1, advances/2);
      }
      else {
        validate(scheme2->isActionRequired(readIterationCheckpoint));
        validateEquals(scheme1->getTimesteps()-1, (advances+1)/2);
      }
    }
    composition.finalize();
    validateEquals(advances, 20);
    validateEquals(scheme1->getTimesteps()-1, 10);
    validateEquals(scheme2->getTimesteps()-1, 10);
  }
  { // Test I(2), E
    preciceDebug("Test I(2), E");
    int numberIterations = 2;
    int maxTimesteps = 10;
    PtrCouplingScheme scheme1(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    numberIterations = 1;
    PtrCouplingScheme scheme2(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    CompositionalCouplingScheme composition;
    composition.addCouplingScheme(scheme1);
    composition.addCouplingScheme(scheme2);
    composition.initialize(0.0, 1);
    int advances = 0;
    while (composition.isCouplingOngoing()){
      composition.advance();
      advances++;
      if (advances%2 == 0){
        validate(scheme1->isActionRequired(writeIterationCheckpoint));
        validateEquals(scheme1->getTimesteps()-1, advances/2);
      }
      else {
        validate(scheme1->isActionRequired(readIterationCheckpoint));
        validateEquals(scheme1->getTimesteps()-1, (advances-1)/2);
      }
    }
    composition.finalize();
    validateEquals(advances, 20);
    validateEquals(scheme1->getTimesteps()-1, 10);
    validateEquals(scheme2->getTimesteps()-1, 10);
  }
  { // Test E, I(3)
    preciceDebug("Test E, I(3)");
    int numberIterations = 1;
    int maxTimesteps = 10;
    PtrCouplingScheme scheme1(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    numberIterations = 3;
    PtrCouplingScheme scheme2(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    CompositionalCouplingScheme composition;
    composition.addCouplingScheme(scheme1);
    composition.addCouplingScheme(scheme2);
    composition.initialize(0.0, 1);
    int advances = 0;
    while (composition.isCouplingOngoing()){
      composition.advance();
      advances++;
      if (advances%3 == 0){
        validate(scheme2->isActionRequired(writeIterationCheckpoint));
        validateEquals(scheme1->getTimesteps()-1, advances/3);
      }
      else {
        validate(scheme2->isActionRequired(readIterationCheckpoint));
        validateEquals(scheme1->getTimesteps()-1, (advances+(3-advances%3))/3);
      }
    }
    composition.finalize();
    validateEquals(advances, 30);
    validateEquals(scheme1->getTimesteps()-1, 10);
    validateEquals(scheme2->getTimesteps()-1, 10);
  }
  { // Test I(3), E
    preciceDebug("Test I(3), E");
    int numberIterations = 3;
    int maxTimesteps = 10;
    PtrCouplingScheme scheme1(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    numberIterations = 1;
    PtrCouplingScheme scheme2(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    CompositionalCouplingScheme composition;
    composition.addCouplingScheme(scheme1);
    composition.addCouplingScheme(scheme2);
    composition.initialize(0.0, 1);
    int advances = 0;
    while (composition.isCouplingOngoing()){
      composition.advance();
      advances++;
      if (advances%3 == 0){
        validate(scheme1->isActionRequired(writeIterationCheckpoint));
        validateEquals(scheme1->getTimesteps()-1, advances/3);
      }
      else {
        validate(scheme1->isActionRequired(readIterationCheckpoint));
        validateEquals(scheme1->getTimesteps()-1, (advances-(advances%3))/3);
      }
    }
    composition.finalize();
    validateEquals(advances, 30);
    validateEquals(scheme1->getTimesteps()-1, 10);
    validateEquals(scheme2->getTimesteps()-1, 10);
  }
  { // Test E, I(2), I(2)
    preciceDebug("Test E, I(2), I(2)");
    int numberIterations = 1;
    int maxTimesteps = 10;
    PtrCouplingScheme scheme1(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    numberIterations = 2;
    PtrCouplingScheme scheme2(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    PtrCouplingScheme scheme3(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    CompositionalCouplingScheme composition;
    composition.addCouplingScheme(scheme1);
    composition.addCouplingScheme(scheme2);
    composition.addCouplingScheme(scheme3);
    composition.initialize(0.0, 1);
    int advances = 0;
    while (composition.isCouplingOngoing()){
      composition.advance();
      advances++;
      if (advances%2 == 0){
        validate(scheme2->isActionRequired(writeIterationCheckpoint));
        validate(scheme3->isActionRequired(writeIterationCheckpoint));
        validateEquals(scheme1->getTimesteps()-1, advances/2);
      }
      else {
        validate(scheme2->isActionRequired(readIterationCheckpoint));
        validate(scheme3->isActionRequired(readIterationCheckpoint));
        validateEquals(scheme1->getTimesteps()-1, (advances+1)/2);
      }
    }
    composition.finalize();
    validateEquals(advances, 20);
    validateEquals(scheme1->getTimesteps()-1, 10);
    validateEquals(scheme2->getTimesteps()-1, 10);
    validateEquals(scheme3->getTimesteps()-1, 10);
  }
  { // Test E, I(2), I(3)
    preciceDebug("Test E, I(2), I(3)");
    int numberIterations = 1;
    int maxTimesteps = 10;
    PtrCouplingScheme scheme1(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    numberIterations = 2;
    PtrCouplingScheme scheme2(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    numberIterations = 3;
    PtrCouplingScheme scheme3(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    CompositionalCouplingScheme composition;
    composition.addCouplingScheme(scheme1);
    composition.addCouplingScheme(scheme2);
    composition.addCouplingScheme(scheme3);
    composition.initialize(0.0, 1);
    int advances = 0;
    while (composition.isCouplingOngoing()){
      composition.advance();
      advances++;
      if (advances%3 == 0){
        validateEquals(scheme1->getTimesteps()-1, advances/3);
        validate(scheme2->isActionRequired(writeIterationCheckpoint));
        validate(scheme3->isActionRequired(writeIterationCheckpoint));
      }
      else if (advances%3 == 1){
        validate(scheme2->isActionRequired(readIterationCheckpoint));
        validate(scheme3->isActionRequired(readIterationCheckpoint));
        validateEquals(scheme1->getTimesteps()-1, (advances+2)/3);
      }
      else if (advances%3 == 2){
        validate(scheme2->isActionRequired(writeIterationCheckpoint));
        validate(scheme3->isActionRequired(readIterationCheckpoint));
        validateEquals(scheme1->getTimesteps()-1, (advances+1)/3);
      }
    }
    composition.finalize();
    validateEquals(advances, 30);
    validateEquals(scheme1->getTimesteps()-1, 10);
    validateEquals(scheme2->getTimesteps()-1, 10);
    validateEquals(scheme3->getTimesteps()-1, 10);
  }
  { // Test I(2), I(2), E
    preciceDebug("Test I(2), I(2), E");
    int numberIterations = 2;
    int maxTimesteps = 10;
    PtrCouplingScheme scheme1(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    PtrCouplingScheme scheme2(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    numberIterations = 1;
    PtrCouplingScheme scheme3(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    CompositionalCouplingScheme composition;
    composition.addCouplingScheme(scheme1);
    composition.addCouplingScheme(scheme2);
    composition.addCouplingScheme(scheme3);
    composition.initialize(0.0, 1);
    int advances = 0;
    while (composition.isCouplingOngoing()){
      composition.advance();
      advances++;
      if (advances%2 == 0){
        validate(scheme1->isActionRequired(writeIterationCheckpoint));
        validate(scheme2->isActionRequired(writeIterationCheckpoint));
        validateEquals(scheme1->getTimesteps()-1, advances/2);
        validateEquals(scheme2->getTimesteps()-1, advances/2);
        validateEquals(scheme3->getTimesteps()-1, advances/2);
      }
      else {
        validate(scheme1->isActionRequired(readIterationCheckpoint));
        validate(scheme2->isActionRequired(readIterationCheckpoint));
        validateEquals(scheme1->getTimesteps()-1, (advances-1)/2);
        validateEquals(scheme2->getTimesteps()-1, (advances-1)/2);
        validateEquals(scheme3->getTimesteps()-1, (advances-1)/2);
      }
    }
    composition.finalize();
    validateEquals(advances, 20);
    validateEquals(scheme1->getTimesteps()-1, 10);
    validateEquals(scheme2->getTimesteps()-1, 10);
    validateEquals(scheme3->getTimesteps()-1, 10);
  }
  { // Test I(2), I(2), E
    preciceDebug("Test I(3), I(2), E");
    int numberIterations = 3;
    int maxTimesteps = 10;
    PtrCouplingScheme scheme1(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    numberIterations = 2;
    PtrCouplingScheme scheme2(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    numberIterations = 1;
    PtrCouplingScheme scheme3(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    CompositionalCouplingScheme composition;
    composition.addCouplingScheme(scheme1);
    composition.addCouplingScheme(scheme2);
    composition.addCouplingScheme(scheme3);
    composition.initialize(0.0, 1);
    int advances = 0;
    while (composition.isCouplingOngoing()){
      composition.advance();
      advances++;
      if (advances%3 == 0){
        validate(scheme1->isActionRequired(writeIterationCheckpoint));
        validate(scheme2->isActionRequired(writeIterationCheckpoint));
        validateEquals(scheme1->getTimesteps()-1, advances/3);
        validateEquals(scheme2->getTimesteps()-1, advances/3);
        validateEquals(scheme3->getTimesteps()-1, advances/3);
      }
      else if (advances%3 == 1){
        validate(scheme1->isActionRequired(readIterationCheckpoint));
        validate(scheme2->isActionRequired(readIterationCheckpoint));
        validateEquals(scheme1->getTimesteps()-1, (advances-1)/3);
        validateEquals(scheme2->getTimesteps()-1, (advances-1)/3);
        validateEquals(scheme3->getTimesteps()-1, (advances-1)/3);
      }
      else if (advances%3 == 2){
        validate(scheme1->isActionRequired(readIterationCheckpoint));
        validate(scheme2->isActionRequired(writeIterationCheckpoint));
        validateEquals(scheme1->getTimesteps()-1, (advances-2)/3);
        validateEquals(scheme2->getTimesteps()-1, (advances+1)/3);
        validateEquals(scheme3->getTimesteps()-1, (advances-2)/3);
      }
    }
    composition.finalize();
    validateEquals(advances, 30);
    validateEquals(scheme1->getTimesteps()-1, 10);
    validateEquals(scheme2->getTimesteps()-1, 10);
    validateEquals(scheme3->getTimesteps()-1, 10);
  }
  {
    preciceDebug("Test I(3), E, I(2)");
    int numberIterations = 3;
    int maxTimesteps = 10;
    PtrCouplingScheme scheme1(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    numberIterations = 1;
    PtrCouplingScheme scheme2(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    numberIterations = 2;
    PtrCouplingScheme scheme3(
      new DummyCouplingScheme(numberIterations, maxTimesteps));
    CompositionalCouplingScheme composition;
    composition.addCouplingScheme(scheme1);
    composition.addCouplingScheme(scheme2);
    composition.addCouplingScheme(scheme3);
    composition.initialize(0.0, 1);
    int advances = 0;
    while (composition.isCouplingOngoing()){
      composition.advance();
      advances++;
      if (advances % 4 >= 3){
        validate(scheme1->isActionRequired(writeIterationCheckpoint));
        validateEquals(scheme1->getTimesteps()-1, (advances-(advances%4)+4)/4);
      }
      else if (advances % 4 != 0){
        validate(scheme1->isActionRequired(readIterationCheckpoint));
      }
      validateEquals(scheme2->getTimesteps()-1, (advances+1)/4);
      validateEquals(scheme2->getTimesteps()-1, (advances+1)/4);
    }
    composition.finalize();
    validateEquals(advances, 40);
    validateEquals(scheme1->getTimesteps()-1, 10);
    validateEquals(scheme2->getTimesteps()-1, 10);
    validateEquals(scheme3->getTimesteps()-1, 10);
  }
}

#ifndef PRECICE_NO_MPI

void CompositionalCouplingSchemeTest:: testExplicitSchemeComposition1()
{
  preciceTrace("testExplicitSchemeComposition1()");
  std::string configPath(_pathToTests + "multi-solver-coupling-1.xml");
  setupAndRunThreeSolverCoupling(configPath);
}

void CompositionalCouplingSchemeTest:: testImplicitSchemeComposition()
{
  preciceTrace("testImplicitSchemeComposition()");
  std::string configPath(_pathToTests + "multi-solver-coupling-2.xml");
  setupAndRunThreeSolverCoupling(configPath);
}

void CompositionalCouplingSchemeTest:: testImplicitExplicitSchemeComposition()
{
  preciceTrace("testImplicitExplicitSchemeComposition()");
  std::string configPath(_pathToTests + "multi-solver-coupling-3.xml");
  setupAndRunThreeSolverCoupling(configPath);
}

void CompositionalCouplingSchemeTest:: testExplicitImplicitSchemeComposition()
{
  preciceTrace("testExplicitImplicitSchemeComposition()");
  std::string configPath(_pathToTests + "multi-solver-coupling-4.xml");
  setupAndRunThreeSolverCoupling(configPath);
}

void CompositionalCouplingSchemeTest:: setupAndRunThreeSolverCoupling
(
  const std::string& configFilename)
{
  preciceTrace1("setupThreeSolverCoupling()", configFilename);
  using namespace mesh;
  utils::Parallel::synchronizeProcesses();
  assertion(utils::Parallel::getCommunicatorSize() > 1);
  mesh::PropertyContainer::resetPropertyIDCounter();

  std::string configurationPath(configFilename);
  std::string nameParticipant0("Participant0");
  std::string nameParticipant1("Participant1");
  std::string nameParticipant2("Participant2");
  std::string localParticipant("");

  utils::XMLTag root = utils::getRootTag();
  PtrDataConfiguration dataConfig(new DataConfiguration(root));
  dataConfig->setDimensions(3);
  PtrMeshConfiguration meshConfig(new MeshConfiguration(root, dataConfig));
  meshConfig->setDimensions(3);
  m2n::M2NConfiguration::SharedPointer m2nConfig(new m2n::M2NConfiguration(root));
  geometry::PtrGeometryConfiguration geoConfig(new geometry::GeometryConfiguration(root, meshConfig));
  geoConfig->setDimensions(3);
  CouplingSchemeConfiguration cplSchemeConfig(root, meshConfig, m2nConfig );

  utils::configure(root, configurationPath);
  meshConfig->setMeshSubIDs();
  m2n::M2N::SharedPointer m2n0 =
      m2nConfig->getM2N(nameParticipant0, nameParticipant1);
  m2n::M2N::SharedPointer m2n1 =
      m2nConfig->getM2N(nameParticipant1, nameParticipant2);

  geoConfig->geometries()[0]->create(*meshConfig->meshes()[0]);

  if (utils::Parallel::getProcessRank() == 0){
    localParticipant = nameParticipant0;
    connect(nameParticipant0, nameParticipant1, localParticipant, m2n0);
  }
  else if (utils::Parallel::getProcessRank() == 1){
    localParticipant = nameParticipant1;
    connect(nameParticipant0, nameParticipant1, localParticipant, m2n0);
    connect(nameParticipant1, nameParticipant2, localParticipant, m2n1);
  }
  else {
    assertion1(utils::Parallel::getProcessRank() == 2,
               utils::Parallel::getProcessRank());
    localParticipant = nameParticipant2;
    connect(nameParticipant1, nameParticipant2, localParticipant, m2n1);
  }

  runThreeSolverCoupling(cplSchemeConfig.getCouplingScheme(localParticipant),
                         localParticipant, meshConfig);
}

void CompositionalCouplingSchemeTest:: runThreeSolverCoupling
(
  PtrCouplingScheme          cplScheme,
  const std::string&         participantName,
  mesh::PtrMeshConfiguration meshConfig )
{
  preciceTrace1("runThreeSolverCoupling()", participantName);

  validateEquals(meshConfig->meshes().size(), 1);
  mesh::PtrMesh mesh = meshConfig->meshes()[0];
  validateEquals(mesh->data().size(), 3);
  validate(mesh->vertices().size() > 0);
  Vector3D valueData1(1.0);
  Vector3D valueData2(1.0);

  std::string readIterationCheckpoint(constants::actionReadIterationCheckpoint());
  std::string writeIterationCheckpoint(constants::actionWriteIterationCheckpoint());

  double computedTime = 0.0;
  int computedTimesteps = 0;

  if (participantName == std::string("Participant0")){
    cplScheme->initialize(0.0, 1);
    validate(not cplScheme->hasDataBeenExchanged());
    validateEquals(cplScheme->isCouplingTimestepComplete(), false);
    validateEquals(cplScheme->isCouplingOngoing(), true);
    while (cplScheme->isCouplingOngoing()){
      validateNumericalEquals(0.1, cplScheme->getNextTimestepMaxLength());
      if (cplScheme->isActionRequired(writeIterationCheckpoint)){
        cplScheme->performedAction(writeIterationCheckpoint);
      }
      cplScheme->addComputedTime(cplScheme->getNextTimestepMaxLength());
      cplScheme->advance();
      if (cplScheme->isActionRequired(readIterationCheckpoint)){
        cplScheme->performedAction(readIterationCheckpoint);
      }
      else {
        validate(cplScheme->isCouplingTimestepComplete());
        computedTime += cplScheme->getNextTimestepMaxLength();
        computedTimesteps++;
      }
      validateNumericalEquals(computedTime, cplScheme->getTime());
      validateEquals(computedTimesteps, cplScheme->getTimesteps()-1);
      validate(cplScheme->hasDataBeenExchanged());
    }
    cplScheme->finalize();
    validateEquals(computedTimesteps, 10);
    validateEquals(cplScheme->isCouplingTimestepComplete(), true);
    validateEquals(cplScheme->isCouplingOngoing(), false);
    validate(cplScheme->getNextTimestepMaxLength() > 0.0); // ??
  }
  else if (participantName == std::string("Participant1")){
    cplScheme->initialize(0.0, 1);
    validate(cplScheme->hasDataBeenExchanged());
    validateEquals(cplScheme->isCouplingTimestepComplete(), false);
    validateEquals(cplScheme->isCouplingOngoing(), true);
    while (cplScheme->isCouplingOngoing()){
      validateNumericalEquals(0.1, cplScheme->getNextTimestepMaxLength());
      if (cplScheme->isActionRequired(writeIterationCheckpoint)){
        cplScheme->performedAction(writeIterationCheckpoint);
      }
      cplScheme->addComputedTime(cplScheme->getNextTimestepMaxLength());
      cplScheme->advance();
      if (cplScheme->isActionRequired(readIterationCheckpoint)){
        cplScheme->performedAction(readIterationCheckpoint);
      }
      else {
        validate(cplScheme->isCouplingTimestepComplete());
        computedTime += cplScheme->getNextTimestepMaxLength();
        computedTimesteps++;
      }
      validateNumericalEquals(computedTime, cplScheme->getTime());
      validateEquals(computedTimesteps, cplScheme->getTimesteps()-1);
      validate(cplScheme->hasDataBeenExchanged());
    }
    cplScheme->finalize();
    validateEquals(computedTimesteps, 10);
    validateEquals(cplScheme->isCouplingTimestepComplete(), true);
    validateEquals(cplScheme->isCouplingOngoing(), false);
    validate(cplScheme->getNextTimestepMaxLength() > 0.0); // ??
  }
  else {
    assertion1(participantName == std::string("Participant2"), participantName);
    cplScheme->initialize(0.0, 1);
    validate(cplScheme->hasDataBeenExchanged());
    validateEquals(cplScheme->isCouplingTimestepComplete(), false);
    validateEquals(cplScheme->isCouplingOngoing(), true);
    while (cplScheme->isCouplingOngoing()){
      validateNumericalEquals(0.1, cplScheme->getNextTimestepMaxLength());
      if (cplScheme->isActionRequired(writeIterationCheckpoint)){
        cplScheme->performedAction(writeIterationCheckpoint);
      }
      cplScheme->addComputedTime(cplScheme->getNextTimestepMaxLength());
      cplScheme->advance();
      if (cplScheme->isActionRequired(readIterationCheckpoint)){
        cplScheme->performedAction(readIterationCheckpoint);
      }
      else {
        validate(cplScheme->isCouplingTimestepComplete());
        computedTime += cplScheme->getNextTimestepMaxLength();
        computedTimesteps++;
      }
      validateNumericalEquals(computedTime, cplScheme->getTime());
      validateEquals(computedTimesteps, cplScheme->getTimesteps()-1);
      validate(cplScheme->hasDataBeenExchanged());
    }
    cplScheme->finalize();
    validateEquals(computedTimesteps, 10);
    validateEquals(cplScheme->isCouplingTimestepComplete(), true);
    validateEquals(cplScheme->isCouplingOngoing(), false);
    validate(cplScheme->getNextTimestepMaxLength() > 0.0); // ??
  }
}

void CompositionalCouplingSchemeTest:: connect
(
  const std::string&     participant0,
  const std::string&     participant1,
  const std::string&     localParticipant,
  m2n::M2N::SharedPointer& communication ) const
{
  preciceTrace3 ( "connect()", participant0, participant1, localParticipant );
  assertion ( communication.use_count() > 0 );
  assertion ( not communication->isConnected() );
  utils::Parallel::initialize ( NULL, NULL, localParticipant );
  if ( participant0 == localParticipant ) {
    communication->requestMasterConnection ( participant1, participant0 );
  }
  else {
    assertion ( participant1 == localParticipant );
    communication->acceptMasterConnection ( participant1, participant0 );
  }
}

#endif // not PRECICE_NO_MPI

}}} // namespace precice, cplscheme, tests
