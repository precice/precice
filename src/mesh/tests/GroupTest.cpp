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

logging::Logger GroupTest:: _log ( "precice::mesh::tests::GroupTest" );

GroupTest:: GroupTest ()
:
   TestCase ("mesh::GroupTest")
{}

void GroupTest:: run()
{
  PRECICE_MASTER_ONLY {
    preciceTrace("run()");
    Group group;
    Vertex vertex0(Eigen::Vector3d::Constant(0.0), 0);
    Vertex vertex1(Eigen::Vector3d::Constant(1.0), 1);
    Vertex vertex2(Eigen::Vector3d::Constant(2.0), 2);
    Vertex vertex3(Eigen::Vector3d::Constant(3.0), 3);
    group.add(vertex0);
    group.add(vertex1);
    group.add(vertex2);
    group.add(vertex3);

    validateEquals(group.size(), 4);

    Eigen::Vector3d coords = Eigen::Vector3d::Zero();
    for (Vertex& v : group.vertices()){
      validate(v.getCoords() == coords);
      coords += Eigen::Vector3d::Constant(1.0);
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

    coords = Eigen::Vector3d::Constant(0.0);
    for (Edge& e : group.edges()){
      validate( e.vertex(0).getCoords() == coords);
      coords += Eigen::Vector3d::Constant(1.0);
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
