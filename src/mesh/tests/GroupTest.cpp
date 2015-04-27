// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "GroupTest.hpp"
#include "../Group.hpp"
#include "../Vertex.hpp"
#include "../Edge.hpp"
#include "../Triangle.hpp"
#include "../Quad.hpp"
#include "utils/Parallel.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Globals.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::mesh::tests::GroupTest)

namespace precice {
namespace mesh {
namespace tests {

tarch::logging::Log GroupTest:: _log ( "precice::mesh::tests::GroupTest" );

GroupTest:: GroupTest ()
:
   TestCase ("mesh::GroupTest")
{}

void GroupTest:: run()
{
  PRECICE_MASTER_ONLY {
    preciceTrace("run()");
    using utils::Vector3D;
    Group group;
    Vertex vertex0(Vector3D(0.0), 0);
    Vertex vertex1(Vector3D(1.0), 1);
    Vertex vertex2(Vector3D(2.0), 2);
    Vertex vertex3(Vector3D(3.0), 3);
    group.add(vertex0);
    group.add(vertex1);
    group.add(vertex2);
    group.add(vertex3);

    validateEquals(group.size(), 4);

    Vector3D coords(0.0);
    for (Vertex& v : group.vertices()){
      validate(tarch::la::equals(v.getCoords(), coords));
      coords += Vector3D(1.0);
    }

    Edge edge0(vertex0, vertex1, 0);
    Edge edge1(vertex1, vertex2, 1);
    Edge edge2(vertex2, vertex3, 2);
    Edge edge3(vertex3, vertex0, 3);
    group.add(edge0);
    group.add(edge1);
    group.add(edge2);
    group.add(edge3);

    validateEquals(group.size(), 8);

    coords = Vector3D(0.0);
    for (Edge& e : group.edges()){
      validate(tarch::la::equals(e.vertex(0).getCoords(), coords));
      coords += Vector3D(1.0);
    }

    Triangle triangle0(edge0, edge1, edge2, 0);
    Triangle triangle1(edge1, edge2, edge3, 1);
    group.add(triangle0);
    group.add(triangle1);

    validateEquals(group.size(), 10);

    Quad quad0(edge0, edge1, edge2, edge3, 0);
    group.add(quad0);

    validateEquals(group.size(), 11);
  }
}


}}} // namespace precice, mesh, tests
