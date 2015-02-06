// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "EdgeTest.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "utils/Parallel.hpp"
#include "utils/Dimensions.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::mesh::tests::EdgeTest)

namespace precice {
namespace mesh {
namespace tests {

tarch::logging::Log EdgeTest:: _log ( "precice::mesh::EdgeTest" );

EdgeTest:: EdgeTest ()
:
   TestCase ("mesh::EdgeTest")
{}

void EdgeTest:: run ()
{
   PRECICE_MASTER_ONLY {
      testMethod ( test );
   }
}

void EdgeTest:: test ()
{
   preciceTrace ( "test()" );

   using utils::Vector3D;
   Vertex v1 ( Vector3D(0.0), 0 );
   Vertex v2 ( Vector3D(1.0), 1 );
//   VertexTuple vertices = { Vertex(Vector(0.0)), Vertex(Vector(1.0)) };
   Edge edge ( v1, v2, 0 );

   Vector3D coords1 = edge.vertex(0).getCoords();
   Vector3D coords2 = edge.vertex(1).getCoords();
   validate ( tarch::la::equals(coords1, Vector3D(0.0)) );
   validate ( tarch::la::equals(coords2, Vector3D(1.0)) );

   Edge edge2 ( v1, v2, 1 );
}

}}} // namespace precice, mesh, tests
