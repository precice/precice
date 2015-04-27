// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "VertexTest.hpp"
#include "mesh/Vertex.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Parallel.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::mesh::tests::VertexTest)

namespace precice {
namespace mesh {
namespace tests {

tarch::logging::Log VertexTest:: _log ( "precice::mesh::VertexTest" );

VertexTest:: VertexTest ()
:
  TestCase ("mesh::VertexTest")
{}

void VertexTest:: run ()
{
  PRECICE_MASTER_ONLY {
    testMethod ( test );
  }
}

void VertexTest:: test ()
{
  preciceTrace ( "test()" );
  using tarch::la::equals;

  Vertex vertex ( utils::Vector3D(1.0), 0 );

  utils::Vector3D coords = vertex.getCoords();
  validate ( tarch::la::equals(coords, utils::Vector3D(1.0)) );

  int id = vertex.getID ();
  validateEquals ( id, 0 );

  utils::Vector3D normal = vertex.getNormal ();
  validate ( tarch::la::equals(normal, utils::Vector3D(0.0)) );

  void* mesh = static_cast<void*> ( vertex.mesh() );
  // Can be replaced by nullptr as soon as we have C++11 available.
  validateEquals ( mesh, static_cast<void*>(NULL) );
}

}}} // namespace precice, mesh, tests
