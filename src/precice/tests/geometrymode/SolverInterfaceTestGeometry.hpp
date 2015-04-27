// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_TESTS_SOLVERINTERFACETESTGEOMETRY_HPP_
#define PRECICE_TESTS_SOLVERINTERFACETESTGEOMETRY_HPP_

#include "precice/SolverInterface.hpp"
#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"
#include "precice/SolverInterface.hpp"

namespace precice {
namespace tests {

/**
 * @brief Provides tests for geometry mode of class SolverInterface.
 */
class SolverInterfaceTestGeometry : public tarch::tests::TestCase
{
public:

  SolverInterfaceTestGeometry();

  virtual ~SolverInterfaceTestGeometry() {}

  /**
   * @brief Retrieves path to test directory.
   */
  virtual void setUp();

  virtual void run();

private:

  static tarch::logging::Log _log;

  std::string _pathToTests;

  int _geoID;

  /**
   * @brief As SolverInterface::configure(), but without changing logging config.
   */
  void configureSolverInterface (
    const std::string& configFilename,
    SolverInterface&   interface );

  void testConfiguration();

  void testManualConfiguration();

  void testSearchQuery();

  void testVoxelQuery();

  void testVoxelQueryMultipleGeometryIDs();

  void testVoxelQueryDFGChannel();

  void testVoxelQueryFSIChannel();

  void testVoxelQueryChannelFour();

  void testVoxelQueryEpsBox();

  void testConservativeStationaryDataMapping();

  /**
   * @brief Tests mapping with radial basis functions.
   */
  void testMappingRBF();

  void testCustomGeometryCreation();

  /**
   * @bried Tests the main functionality that is necessary to perform a Pinelli-type Direct Forcing method
   */
  void testPinelli();

  void testSetExportLocation();

  void testSpacetrees();

  void testBug();

  void testBug2();

  /**
   * @brief Geometry query bug with Peano and static octree involved.
   */
  void testBug3();

  /**
   * @brief Geometry query bug with Peano and dynamic octree involved.
   */
  void testBug4();

  /**
   * @brief Geometry query bug with Peano and dynamic-octree involved.
   */
  void testBug5();

  /**
   * @brief Geometry query bug with Peano and dynamic-octree involved.
   */
  //void testBug6(); No bug, actually

  /**
   * @brief Tests if configured data actions are issued and performed correctly.
   */
  void testDataActions();

  /**
   * @brief Tests updating a spacetree due to mesh modifications.
   */
  void testUpdateSpacetree();

  /**
   * @brief Tests spacetree with multiple meshes.
   */
  void testMultipleMeshSpacetree();
};

}} // namespace precice, tests

#endif /* PRECICE_TESTS_SOLVERINTERFACETESTGEOMETRY_HPP_ */
