#include "DynamicOctreeTest.hpp"
#include "spacetree/ExportSpacetree.hpp"
#include "spacetree/config/SpacetreeConfiguration.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "io/ExportVTK.hpp"
#include "query/ExportVTKNeighbors.hpp"
#include "query/ExportVTKVoxelQueries.hpp"
#include "query/FindVoxelContent.hpp"
#include "geometry/Sphere.hpp"
#include "geometry/Cuboid.hpp"
#include "precice/Constants.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"
#include <iostream>

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::spacetree::tests::DynamicOctreeTest)

namespace precice {
namespace spacetree {
namespace tests {

logging::Logger DynamicOctreeTest::_log ( "precice::spacetree::tests::DynamicOctreeTest" );

DynamicOctreeTest:: DynamicOctreeTest()
:
  TestCase ( "spacetree::DynamicOctreeTest" )
{}

void DynamicOctreeTest:: run()
{
  TRACE();
  int sizeRanks = utils::Parallel::getCommunicatorSize();
  bool evenTasks = false;
  bool oddTasks = false;
  if (sizeRanks > 1){
    if (utils::Parallel::getProcessRank() == 0){
      evenTasks = true;
    }
    else if (utils::Parallel::getProcessRank() == 1){
      oddTasks = true;
    }
  }
  else {
    evenTasks = true;
    oddTasks = true;
  }
  if (evenTasks || oddTasks){
    std::string testName("DynamicOctreeTest");
    DynamicOctreeFactory factory;
    SpacetreeTestScenarios testScenarios(testName, factory);

    if (evenTasks){
      DEBUG ( "test search position" );
      testScenarios.testSearchPosition();
    }

    if (oddTasks){
      DEBUG ( "test search distance" );
      testScenarios.testSearchDistance();
    }

    if (evenTasks){
      DEBUG ( "test neighbor search" );
      testScenarios.testNeighborSearch();
    }

    if (oddTasks){
      DEBUG ( "test search content vertices" );
      testScenarios.testSearchContentVertices();
    }

    if (evenTasks){
      DEBUG ( "test search content edges" );
      testScenarios.testSearchContentEdges();
    }

    if (oddTasks){
      DEBUG ( "test search content triangles" );
      testScenarios.testSearchContentTriangles();
    }

    if (evenTasks){
      DEBUG ( "test search voxel position" );
      testScenarios.testVoxelPosition();
    }

    if (oddTasks){
      DEBUG ( "test splitting voxels" );
      testScenarios.testSplittingVoxels();
    }
    validate(testScenarios.getNumberOfErrors() == 0);
  }
}

}}} // namespace precice, spacetree, tests



