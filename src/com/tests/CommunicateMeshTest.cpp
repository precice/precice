#include "CommunicateMeshTest.hpp"
#include "com/CommunicateMesh.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/PropertyContainer.hpp"
#include "utils/Parallel.hpp"
#include "utils/Helpers.hpp"
#include "math/math.hpp"


#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::com::tests::CommunicateMeshTest)

namespace precice {
namespace com {
namespace tests {

logging::Logger CommunicateMeshTest::
  _log ( "precice::com::tests::CommunicateMeshTest" );

CommunicateMeshTest:: CommunicateMeshTest ()
:
  TestCase ( "precice::com::tests::CommunicateMeshTest" )
{}

void CommunicateMeshTest:: run ()
{
# ifndef PRECICE_NO_MPI
  if ( utils::Parallel::getCommunicatorSize() >= 2 ) {
    testMethod ( testTwoSolvers );
  }
# endif
}

#ifndef PRECICE_NO_MPI

void CommunicateMeshTest:: testTwoSolvers ()
{
  TRACE();
  utils::Parallel::synchronizeProcesses ();
  assertion ( utils::Parallel::getCommunicatorSize() > 1 );
  mesh::PropertyContainer::resetPropertyIDCounter ();

  std::string participant0 ( "rank0" );
  std::string participant1 ( "rank1" );

  for ( int dim=2; dim <= 3; dim++ ){
    // Build mesh to communicate for rank0
    mesh::Mesh mesh ( "MyMesh", dim, false );
    if ( utils::Parallel::getProcessRank() == 0 ){
      mesh::Vertex& v0 = mesh.createVertex ( Eigen::VectorXd::Constant(dim, 0) );
      mesh::Vertex& v1 = mesh.createVertex ( Eigen::VectorXd::Constant(dim, 1) );
      mesh::Vertex& v2 = mesh.createVertex ( Eigen::VectorXd::Constant(dim, 2) );

      mesh.createEdge ( v0, v1 );
      mesh.createEdge ( v1, v2 );
      mesh.createEdge ( v2, v0 );
    }

    // Create mesh communicator
    std::vector<int> involvedRanks = {0, 1};
    MPI_Comm comm = utils::Parallel::getRestrictedCommunicator ( involvedRanks );
    if ( utils::Parallel::getProcessRank() < 2 ) {
      utils::Parallel::setGlobalCommunicator ( comm );
      validateEquals ( utils::Parallel::getCommunicatorSize(), 2 );
      com::Communication::SharedPointer com ( new com::MPIDirectCommunication() );
      CommunicateMesh comMesh ( com );

      if ( utils::Parallel::getProcessRank() == 0 ) {
        utils::Parallel::splitCommunicator(participant0 );
        com->acceptConnection ( participant0, participant1, 0, 1 );
        comMesh.sendMesh ( mesh, 0 );
        validateEquals ( mesh.vertices().size(), 3 );
        validateEquals ( mesh.edges().size(), 3 );
        validate ( math::equals(mesh.vertices()[0].getCoords(), Eigen::VectorXd::Constant(dim, 0) ));
        validate ( math::equals(mesh.vertices()[1].getCoords(), Eigen::VectorXd::Constant(dim, 1) ));
        validate ( math::equals(mesh.vertices()[2].getCoords(), Eigen::VectorXd::Constant(dim, 2) ));
      }
      else if ( utils::Parallel::getProcessRank() == 1 ) {
        mesh.createVertex ( Eigen::VectorXd::Constant(dim, 9) ); // new version receiveMesh can also deal with delta meshes
        utils::Parallel::splitCommunicator(participant1 );
        com->requestConnection ( participant0, participant1, 0, 1 );
        comMesh.receiveMesh ( mesh, 0 );
        validateEquals ( mesh.vertices().size(), 4 );
        validateEquals ( mesh.edges().size(), 3 );
        validate ( math::equals(mesh.vertices()[0].getCoords(), Eigen::VectorXd::Constant(dim, 9) ));
        validate ( math::equals(mesh.vertices()[1].getCoords(), Eigen::VectorXd::Constant(dim, 0) ));
        validate ( math::equals(mesh.vertices()[2].getCoords(), Eigen::VectorXd::Constant(dim, 1) ));
        validate ( math::equals(mesh.vertices()[3].getCoords(), Eigen::VectorXd::Constant(dim, 2) ));

      }
      com->closeConnection ();

      validate ( math::equals(mesh.edges()[0].vertex(0).getCoords(), Eigen::VectorXd::Constant(dim, 0) ));
      validate ( math::equals(mesh.edges()[0].vertex(1).getCoords(), Eigen::VectorXd::Constant(dim, 1) ));
      validate ( math::equals(mesh.edges()[1].vertex(0).getCoords(), Eigen::VectorXd::Constant(dim, 1) ));
      validate ( math::equals(mesh.edges()[1].vertex(1).getCoords(), Eigen::VectorXd::Constant(dim, 2) ));
      validate ( math::equals(mesh.edges()[2].vertex(0).getCoords(), Eigen::VectorXd::Constant(dim, 2) ));
      validate ( math::equals(mesh.edges()[2].vertex(1).getCoords(), Eigen::VectorXd::Constant(dim, 0) ));
      utils::Parallel::clearGroups();
      utils::Parallel::setGlobalCommunicator(utils::Parallel::getCommunicatorWorld());
    }
  }
}

#endif // not PRECICE_NO_MPI

}}} // namespace precice, com, tests
