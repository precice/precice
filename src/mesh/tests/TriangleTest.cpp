// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "TriangleTest.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "utils/Parallel.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Globals.hpp"
#include "boost/array.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::mesh::tests::TriangleTest)

namespace precice {
namespace mesh {
namespace tests {

tarch::logging::Log TriangleTest:: _log ("precice::mesh::TriangleTest");

TriangleTest:: TriangleTest ()
:
 TestCase ("mesh::TriangleTest")
{}

void TriangleTest:: run ()
{
  PRECICE_MASTER_ONLY {
    testMethod ( test );
  }
}

void TriangleTest:: test ()
{
  preciceTrace ( "test" );
  using utils::Vector3D;
  Vector3D coords1 ( 0.0, 0.0, 0.0 );
  Vector3D coords2 ( 1.0, 0.0, 0.0 );
  Vector3D coords3 ( 1.0, 1.0, 0.0 );
  {
    Vertex v1 ( coords1, 0 );
    Vertex v2 ( coords2, 1 );
    Vertex v3 ( coords3, 2 );

    Edge e1 ( v1, v2, 0 );
    Edge e2 ( v3, v2, 1 );
    Edge e3 ( v3, v1, 2 );

    Triangle triangle ( e1, e2, e3, 0 );

    Vertex& v1ref = triangle.vertex ( 0 );
    validateEquals ( v1ref.getID(), v1.getID() );

    Vertex& v2ref = triangle.vertex ( 1 );
    validateEquals ( v2ref.getID(), v2.getID() );

    Vertex& v3ref = triangle.vertex ( 2 );
    validateEquals ( v3ref.getID(), v3.getID() );

    Edge& e1ref = triangle.edge ( 0 );
    validateEquals ( e1ref.getID(), e1.getID() );

    Edge& e2ref = triangle.edge ( 1 );
    validateEquals ( e2ref.getID(), e2.getID() );

    Edge& e3ref = triangle.edge ( 2 );
    validateEquals ( e3ref.getID(), e3.getID() );

    int id = triangle.getID ();
    validateEquals ( id, 0 );
  }
  {
    Vertex v1 ( coords1, 0 );
    Vertex v2 ( coords2, 1 );
    Vertex v3 ( coords3, 2 );

    Edge e1 ( v1, v2, 0 );
    Edge e2 ( v3, v2, 1 );
    Edge e3 ( v1, v3, 2 );

    Triangle triangle ( e1, e2, e3, 0 );

    Vertex& v1ref = triangle.vertex ( 0 );
    validateEquals ( v1ref.getID(), v1.getID() );

    Vertex& v2ref = triangle.vertex ( 1 );
    validateEquals ( v2ref.getID(), v2.getID() );

    Vertex& v3ref = triangle.vertex ( 2 );
    validateEquals ( v3ref.getID(), v3.getID() );

    Edge& e1ref = triangle.edge ( 0 );
    validateEquals ( e1ref.getID(), e1.getID() );

    Edge& e2ref = triangle.edge ( 1 );
    validateEquals ( e2ref.getID(), e2.getID() );

    Edge& e3ref = triangle.edge ( 2 );
    validateEquals ( e3ref.getID(), e3.getID() );

    int id = triangle.getID ();
    validateEquals ( id, 0 );
  }
  {
    Vertex v1 ( coords1, 0 );
    Vertex v2 ( coords2, 1 );
    Vertex v3 ( coords3, 2 );

    Edge e1 ( v1, v2, 0 );
    Edge e2 ( v3, v2, 1 );
    Edge e3 ( v3, v1, 2 );

    Triangle triangle ( e1, e3, e2, 0 );

    Vertex& v1ref = triangle.vertex ( 0 );
    validateEquals ( v1ref.getID(), v2.getID() );

    Vertex& v2ref = triangle.vertex ( 1 );
    validateEquals ( v2ref.getID(), v1.getID() );

    Vertex& v3ref = triangle.vertex ( 2 );
    validateEquals ( v3ref.getID(), v3.getID() );

    Edge& e1ref = triangle.edge ( 0 );
    validateEquals ( e1ref.getID(), e1.getID() );

    Edge& e2ref = triangle.edge ( 1 );
    validateEquals ( e2ref.getID(), e3.getID() );

    Edge& e3ref = triangle.edge ( 2 );
    validateEquals ( e3ref.getID(), e2.getID() );

    int id = triangle.getID ();
    validateEquals ( id, 0 );
  }
}

}}} // namespace precice, mesh, tests
