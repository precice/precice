#ifndef PRECICE_SPACETREE_DYNAMICOCTREETEST_HPP_
#define PRECICE_SPACETREE_DYNAMICOCTREETEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"
#include "utils/Dimensions.hpp"
#include "spacetree/DynamicOctree.hpp"
#include "spacetree/tests/SpacetreeTestScenarios.hpp"

namespace precice {
namespace spacetree {
namespace tests {

/**
 * @brief Provides test cases for class RegularSpacetree.
 */
class DynamicOctreeTest : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor.
   */
  DynamicOctreeTest();

  /**
   * @brief Empty.
   */
  virtual void setUp() {}

  /**
   * @brief Runs all tests.
   */
	virtual void run();

private:

  struct DynamicOctreeFactory : public SpacetreeTestScenarios::SpacetreeFactory
  {
    virtual PtrSpacetree createSpacetree (
      const utils::DynVector& offset,
      const utils::DynVector& halflengths,
      double                  refinementLimit )
    {
      return PtrSpacetree(new DynamicOctree(offset, halflengths[0], refinementLimit));
    }
  };

	// @brief Logging device.
	static tarch::logging::Log _log;

//	/**
//	 * Tests method searchPosition().
//	 */
//	void testSearchPosition();
//
//	void testSearchDistance();
//
//	/**
//	 * @brief Tests finding the closest neighboring point by exporting the
//	 * results given by the accept() and the acceptSearch() methods.
//	 */
//	void testNeighborSearch();
//
//	/**
//	 * @brief Triggers performTestSearchContentVertices() with various parameters.
//	 */
//	void testSearchContentVertices();
//
//	 /**
//	   * @brief Tests finding vertices with the searchContent() method.
//	   */
//	void performTestSearchContentVertices (
//	  int                     dimension,
//	  bool                    positive,
//	  const utils::DynVector& offset);
//
//  /**
//   * @brief Triggers performTestSearchContentEdges() with various parameters.
//   */
//	void testSearchContentEdges();
//
//  /**
//   * @brief Tests finding edges with the searchContent() method.
//   */
//	void performTestSearchContentEdges (
//	  int                     dimension,
//    bool                    positive,
//    const utils::DynVector& offset );
//
//  /**
//   * @brief Triggers performTestSearchContentTriangles() with various parameters.
//   */
//	void testSearchContentTriangles();
//
//  /**
//   * @brief Tests finding triangles with the searchContent() method.
//   */
//	void performTestSearchContentTriangles (
//    int  dimension,
//    int  secondDimension,
//    bool positive );
//
//	/**
//	 * @brief Tests the correctness of some voxel positions obtained by using the
//	 * acceptVoxel() method. The geometry is  a cuboid.
//	 */
//	void testVoxelPosition();
//
//
//	/**
//	 * @brief Tests the correctness of some voxel positions obtained by splitting
//	 * the voxels by callling the acceptVoxel() method.
//	 * The geometry is  a cuboid.
//	 */
//	void testSplittingVoxels();
//
//	/**
//	 * @brief Tests reading xml files with configurations for spacetrees.
//	 */
//	void testConfiguration();
};

}}} // namespace precice, spacetree, tests

#endif /* PRECICE_SPACETREE_DYNAMICOCTREETEST_HPP_ */
