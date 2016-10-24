#include "StaticOctreeTest.hpp"
#include "utils/Helpers.hpp"
#include "utils/Parallel.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::spacetree::tests::StaticOctreeTest)

namespace precice {
namespace spacetree {
namespace tests {

logging::Logger StaticOctreeTest::_log("precice::spacetree::tests::StaticOctreeTest");

void StaticOctreeTest:: run()
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
    std::string testName("StaticOctreeTest");
    StaticOctreeFactory factory;
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
