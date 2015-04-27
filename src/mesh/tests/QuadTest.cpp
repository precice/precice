// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "QuadTest.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Quad.hpp"
#include "utils/Parallel.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Globals.hpp"
#include "boost/array.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::mesh::tests::QuadTest)

namespace precice {
namespace mesh {
namespace tests {

tarch::logging::Log QuadTest:: _log ("precice::mesh::QuadTest");

QuadTest:: QuadTest ()
:
 TestCase ("mesh::QuadTest")
{}

void QuadTest:: run ()
{
  PRECICE_MASTER_ONLY {
    testMethod(test);
  }
}

void QuadTest:: test()
{
  preciceTrace("test");
  using utils::Vector3D;
  Vector3D coords0(0.0, 0.0, 0.0);
  Vector3D coords1(1.0, 0.0, 0.0);
  Vector3D coords2(1.0, 1.0, 0.0);
  Vector3D coords3(0.0, 1.0, 0.0);
  {
    Vertex v0(coords0, 0);
    Vertex v1(coords1, 1);
    Vertex v2(coords2, 2);
    Vertex v3(coords3, 3);

    Edge e0(v0, v1, 0);
    Edge e1(v1, v2, 1);
    Edge e2(v2, v3, 2);
    Edge e3(v3, v0, 3);

    Quad quad(e0, e1, e2, e3, 0);

    Vertex& v0ref = quad.vertex(0);
    validateEquals(v0ref.getID(), v0.getID());
    Vertex& v1ref = quad.vertex(1);
    validateEquals(v1ref.getID(), v1.getID());
    Vertex& v2ref = quad.vertex(2);
    validateEquals(v2ref.getID(), v2.getID());
    Vertex& v3ref = quad.vertex(3);
    validateEquals(v3ref.getID(), v3.getID());

    Edge& e0ref = quad.edge(0);
    validateEquals(e0ref.getID(), e0.getID());
    Edge& e1ref = quad.edge(1);
    validateEquals (e1ref.getID(), e1.getID());
    Edge& e2ref = quad.edge(2);
    validateEquals(e2ref.getID(), e2.getID());
    Edge& e3ref = quad.edge(3);
    validateEquals(e3ref.getID(), e3.getID());


    int id = quad.getID();
    validateEquals(id, 0);
  }
  {
    Vertex v0(coords0, 0);
    Vertex v1(coords1, 1);
    Vertex v2(coords2, 2);
    Vertex v3(coords3, 3);

    Edge e0(v0, v1, 0);
    Edge e1(v1, v2, 1);
    Edge e2(v3, v2, 2);
    Edge e3(v0, v3, 3);

    Quad quad(e0, e1, e2, e3, 0);

    Vertex& v0ref = quad.vertex(0);
    validateEquals(v0ref.getID(), v0.getID());
    Vertex& v1ref = quad.vertex(1);
    validateEquals(v1ref.getID(), v1.getID());
    Vertex& v2ref = quad.vertex(2);
    validateEquals(v2ref.getID(), v2.getID());
    Vertex& v3ref = quad.vertex(3);
    validateEquals(v3ref.getID(), v3.getID());

    Edge& e0ref = quad.edge(0);
    validateEquals(e0ref.getID(), e0.getID());
    Edge& e1ref = quad.edge(1);
    validateEquals(e1ref.getID(), e1.getID());
    Edge& e2ref = quad.edge(2);
    validateEquals(e2ref.getID(), e2.getID());
    Edge& e3ref = quad.edge(3);
    validateEquals(e3ref.getID(), e3.getID());

    int id = quad.getID();
    validateEquals(id, 0);
  }
  {
    Vertex v0(coords0, 0);
    Vertex v1(coords1, 1);
    Vertex v2(coords2, 2);
    Vertex v3(coords3, 3);

    Edge e0(v0, v1, 0);
    Edge e1(v2, v1, 1);
    Edge e2(v2, v3, 2);
    Edge e3(v0, v3, 3);

    Quad quad(e0, e1, e2, e3, 0);

    Vertex& v0ref = quad.vertex(0);
    validateEquals(v0ref.getID(), v0.getID());
    Vertex& v1ref = quad.vertex(1);
    validateEquals(v1ref.getID(), v1.getID());
    Vertex& v2ref = quad.vertex(2);
    validateEquals(v2ref.getID(), v2.getID());
    Vertex& v3ref = quad.vertex(3);
    validateEquals(v3ref.getID(), v3.getID());

    Edge& e0ref = quad.edge(0);
    validateEquals(e0ref.getID(), e0.getID());
    Edge& e1ref = quad.edge(1);
    validateEquals(e1ref.getID(), e1.getID());
    Edge& e2ref = quad.edge(2);
    validateEquals(e2ref.getID(), e2.getID());
    Edge& e3ref = quad.edge(3);
    validateEquals(e3ref.getID(), e3.getID());

    int id = quad.getID();
    validateEquals(id, 0);
  }
}

}}} // namespace precice, mesh, tests
