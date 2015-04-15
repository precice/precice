// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "CommunicateMeshTest.hpp"
#include "com/CommunicateMesh.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/PropertyContainer.hpp"
#include "utils/Parallel.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Helpers.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::com::tests::CommunicateMeshTest)

namespace precice {
namespace com {
namespace tests {

tarch::logging::Log CommunicateMeshTest::
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
  preciceTrace ( "testTwoSolvers" );
  using utils::DynVector;
  utils::Parallel::synchronizeProcesses ();
  assertion ( utils::Parallel::getCommunicatorSize() > 1 );
  mesh::PropertyContainer::resetPropertyIDCounter ();

  std::string participant0 ( "rank0" );
  std::string participant1 ( "rank1" );

  for ( int dim=2; dim <= 3; dim++ ){
    // Build mesh to communicate for rank0
    mesh::Mesh mesh ( "MyMesh", dim, false );
    if ( utils::Parallel::getProcessRank() == 0 ){
      mesh::Vertex& v0 = mesh.createVertex ( DynVector(dim,0.0) );
      mesh::Vertex& v1 = mesh.createVertex ( DynVector(dim,1.0) );
      mesh::Vertex& v2 = mesh.createVertex ( DynVector(dim,2.0) );

      mesh.createEdge ( v0, v1 );
      mesh.createEdge ( v1, v2 );
      mesh.createEdge ( v2, v0 );
    }

    // Create mesh communicator
    std::vector<int> involvedRanks;
    involvedRanks += 0, 1;
    MPI_Comm comm = utils::Parallel::getRestrictedCommunicator ( involvedRanks );
    if ( utils::Parallel::getProcessRank() < 2 ) {
      utils::Parallel::setGlobalCommunicator ( comm );
      validateEquals ( utils::Parallel::getCommunicatorSize(), 2 );
      com::Communication::SharedPointer com ( new com::MPIDirectCommunication() );
      CommunicateMesh comMesh ( com );

      if ( utils::Parallel::getProcessRank() == 0 ) {
        utils::Parallel::initialize ( NULL, NULL, participant0 );
        com->acceptConnection ( participant0, participant1, 0, 1 );
        comMesh.sendMesh ( mesh, 0 );
        validateEquals ( mesh.vertices().size(), 3 );
        validateEquals ( mesh.edges().size(), 3 );
        validate ( equals(mesh.vertices()[0].getCoords(), DynVector(dim,0.0)) );
        validate ( equals(mesh.vertices()[1].getCoords(), DynVector(dim,1.0)) );
        validate ( equals(mesh.vertices()[2].getCoords(), DynVector(dim,2.0)) );
      }
      else if ( utils::Parallel::getProcessRank() == 1 ) {
        mesh.createVertex ( DynVector(dim,9.0) ); // new version receiveMesh can also deal with delta meshes
        utils::Parallel::initialize ( NULL, NULL, participant1 );
        com->requestConnection ( participant0, participant1, 0, 1 );
        comMesh.receiveMesh ( mesh, 0 );
        validateEquals ( mesh.vertices().size(), 4 );
        validateEquals ( mesh.edges().size(), 3 );
        validate ( equals(mesh.vertices()[0].getCoords(), DynVector(dim,9.0)) );
        validate ( equals(mesh.vertices()[1].getCoords(), DynVector(dim,0.0)) );
        validate ( equals(mesh.vertices()[2].getCoords(), DynVector(dim,1.0)) );
        validate ( equals(mesh.vertices()[3].getCoords(), DynVector(dim,2.0)) );

      }
      com->closeConnection ();

      using tarch::la::equals;

      validate ( equals(mesh.edges()[0].vertex(0).getCoords(), DynVector(dim,0.0)) );
      validate ( equals(mesh.edges()[0].vertex(1).getCoords(), DynVector(dim,1.0)) );
      validate ( equals(mesh.edges()[1].vertex(0).getCoords(), DynVector(dim,1.0)) );
      validate ( equals(mesh.edges()[1].vertex(1).getCoords(), DynVector(dim,2.0)) );
      validate ( equals(mesh.edges()[2].vertex(0).getCoords(), DynVector(dim,2.0)) );
      validate ( equals(mesh.edges()[2].vertex(1).getCoords(), DynVector(dim,0.0)) );
      utils::Parallel::setGlobalCommunicator(utils::Parallel::getCommunicatorWorld());
    }
  }
}

#endif // not PRECICE_NO_MPI

}}} // namespace precice, com, tests
