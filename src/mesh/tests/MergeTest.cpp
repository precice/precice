// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "MergeTest.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Quad.hpp"
#include "mesh/Group.hpp"
#include "mesh/Merge.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::mesh::tests::MergeTest)

namespace precice {
namespace mesh {
namespace tests {

tarch::logging::Log MergeTest:: _log ( "precice::mesh::tests::MergeTest" );

MergeTest:: MergeTest()
:
   TestCase ("mesh::MergeTest")
{}

void MergeTest:: run ()
{
  PRECICE_MASTER_ONLY {
    preciceTrace("run()");
    using utils::Vector3D;
    // Create visitables
    int dim = 3;
    Mesh mesh("MyMesh", dim, false);
    Vertex& v0 = mesh.createVertex(utils::Vector3D(0.0));
    Vertex& v1 = mesh.createVertex(utils::Vector3D(0.0));
    Vertex& v2 = mesh.createVertex(utils::Vector3D(0.0));
    Vertex& v3 = mesh.createVertex(utils::Vector3D(0.0));
    Edge& e0 = mesh.createEdge(v0, v1);
    Edge& e1 = mesh.createEdge(v1, v2);
    Edge& e2 = mesh.createEdge(v2, v0);
    Edge& e3 = mesh.createEdge(v2, v3);
    Edge& e4 = mesh.createEdge(v3, v0);
    Triangle& t1 = mesh.createTriangle(e0, e1, e2);
    Triangle& t2 = mesh.createTriangle(e2, e1, e0);
    Quad& q0 = mesh.createQuad(e0, e1, e3, e4);

    Group group;
    group.add(v1);
    group.add(v2);
    group.add(v3);
    group.add(v2);
    group.add(v1);
    group.add(v3);
    group.add(v2);

    group.add(e1);
    group.add(e1);
    group.add(e2);
    group.add(e3);
    group.add(e3);
    group.add(e4);

    group.add(t1);
    group.add(t1);
    group.add(t2);
    group.add(t2);

    group.add(q0);
    group.add(q0);

    Merge merge;
    merge(group);
    validateEquals(merge.content().vertices().size(), 3);
    validateEquals(merge.content().edges().size(), 4);
    validateEquals(merge.content().triangles().size(), 2);
    validateEquals(merge.content().quads().size(), 1);
    validateEquals(merge.content().size(), 10);
    merge(group); // Shouldn't change anything
    validateEquals(merge.content().vertices().size(), 3);
    validateEquals(merge.content().edges().size(), 4);
    validateEquals(merge.content().triangles().size(), 2);
    validateEquals(merge.content().quads().size(), 1);
    validateEquals(merge.content().size(), 10);
  }
}

}}} // namespace precice, mesh, tests
