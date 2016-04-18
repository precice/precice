// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_SPACETREE_TESTS_SPACETREETESTSCENARIOS_HPP_
#define PRECICE_SPACETREE_TESTS_SPACETREETESTSCENARIOS_HPP_

#include "spacetree/Spacetree.hpp"
#include "utils/Dimensions.hpp"
#include "tarch/tests/TestCase.h"
#include <string>

namespace precice {
namespace spacetree {
namespace tests {

class SpacetreeTestScenarios : public tarch::tests::TestCase
{
public:

  struct SpacetreeFactory
  {
    virtual PtrSpacetree createSpacetree (
      const utils::DynVector& offset,
      const utils::DynVector& halflengths,
      double                  refinementLimit ) =0;
  };

  SpacetreeTestScenarios (
    const std::string& testName,
    SpacetreeFactory&  factory );

  virtual void setUp() {}

  virtual void run() {}

  /**
   * Tests method searchPosition().
   */
  void testSearchPosition();

  void testSearchDistance();

  /**
   * @brief Tests finding the closest neighboring point by exporting the
   * results given by the accept() and the acceptSearch() methods.
   */
  void testNeighborSearch();

  /**
   * @brief Triggers performTestSearchContentVertices() with various parameters.
   */
  void testSearchContentVertices();

  /**
   * @brief Triggers performTestSearchContentEdges() with various parameters.
   */
  void testSearchContentEdges();

  /**
   * @brief Triggers performTestSearchContentTriangles() with various parameters.
   */
  void testSearchContentTriangles();

  /**
   * @brief Tests the correctness of some voxel positions obtained by using the
   * acceptVoxel() method. The geometry is  a cuboid.
   */
  void testVoxelPosition();

  /**
   * @brief Tests the correctness of some voxel positions obtained by splitting
   * the voxels by callling the acceptVoxel() method.
   * The geometry is  a cuboid.
   */
  void testSplittingVoxels();

private:

  static logging::Logger _log;

  std::string _testName;

  SpacetreeFactory& _factory;

  /**
    * @brief Tests finding vertices with the searchContent() method.
    */
 void performTestSearchContentVertices (
   int                     dimension,
   bool                    positive,
   const utils::DynVector& offset);

 /**
  * @brief Tests finding edges with the searchContent() method.
  */
 void performTestSearchContentEdges (
   int                     dimension,
   bool                    positive,
   const utils::DynVector& offset );

 /**
  * @brief Tests finding triangles with the searchContent() method.
  */
 void performTestSearchContentTriangles (
   int  dimension,
   int  secondDimension,
   bool positive );
};

}}} // namespace precice, spacetree, tests

#endif /* PRECICE_SPACETREE_TESTS_SPACETREETESTSCENARIOS_HPP_ */
