// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "FindClosestVertexVisitorTest.hpp"
#include "query/FindClosestVertex.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "utils/Parallel.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::query::tests::FindClosestVertexVisitorTest)

namespace precice {
namespace query {
namespace tests {

tarch::logging::Log FindClosestVertexVisitorTest::
   _log ("precice::query::FindClosestVertexVisitorTest");

FindClosestVertexVisitorTest:: FindClosestVertexVisitorTest ()
:
  tarch::tests::TestCase ("query::FindClosestVertexVisitorTest")
{}

void FindClosestVertexVisitorTest:: run ()
{
  PRECICE_MASTER_ONLY {
    preciceTrace ( "run()" );
    using utils::Vector2D;
    mesh::Mesh mesh ( "Mesh", 2, false );
    mesh.createVertex ( Vector2D(0.0, 0.0) );
    mesh.createVertex ( Vector2D(0.0, 5.0) );
    FindClosestVertex find ( Vector2D(1.0, 0.0) );
    bool found = find ( mesh );
    validate ( found );
    mesh::Vertex& closestVertex = find.getClosestVertex();
    validate ( tarch::la::equals(closestVertex.getCoords(), Vector2D(0.0,0.0)) );
    double distance = find.getEuclidianDistance ();
    validateEquals (distance, 1.0);
  }
}

}}} // namespace precice, query, tests
