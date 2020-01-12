#include "precice/MeshHandle.hpp"
#include "testing/Testing.hpp"
#include "mesh/Mesh.hpp"
#include <Eigen/Core>
#include <vector>
#include <algorithm>

using namespace precice;

struct MeshFixture {
    MeshFixture() : mesh("MyMesh", 3, false, testing::nextMeshID()), handle(mesh) {
        auto & v1 = mesh.createVertex(Eigen::Vector3d(0, 2, 0));
        auto & v2 = mesh.createVertex(Eigen::Vector3d(2, 1, 0));
        auto & v3 = mesh.createVertex(Eigen::Vector3d(3, 0, 0));
        auto & v4 = mesh.createVertex(Eigen::Vector3d(1, 0, 0));
        // Quad Borders
        auto & e1 = mesh.createEdge(v1, v2);
        auto & e2 = mesh.createEdge(v2, v3);
        auto & e3 = mesh.createEdge(v3, v4);
        auto & e4 = mesh.createEdge(v4, v1);
        // Diagonal
        auto & e5 = mesh.createEdge(v2, v4);
        // Triangles
        mesh.createTriangle(e1, e5, e4);
        mesh.createTriangle(e2, e3, e5);
        // Quad
        mesh.createQuad(e1, e2, e3, e4);
        // Check the Mesh
        BOOST_TEST(mesh.vertices().size() == 4);
        BOOST_TEST(mesh.edges().size() == 5);
        BOOST_TEST(mesh.triangles().size() == 2);
        BOOST_TEST(mesh.quads().size() == 1);
    }

    mesh::Mesh mesh;
    MeshHandle handle;
    const int vertex_cnt = 4;
    const int edge_cnt = 5;
    const int triangle_cnt = 2;
    const int quad_cnt = 1;
    const int primitive_cnt = 4+5+2+1;
};


BOOST_FIXTURE_TEST_SUITE(PreciceTests, MeshFixture)

BOOST_AUTO_TEST_CASE(VertexHandle)
{
    const auto& vertices = handle.vertices();
    BOOST_TEST(std::distance(vertices.begin(), vertices.end()) == vertex_cnt);
    BOOST_TEST(vertices.size() == vertex_cnt);
    std::set<int> vs;
    for(const auto & v : vertices) {
        vs.insert(v.vertexID());
    }
    BOOST_TEST(vs.size() == vertex_cnt);
}

BOOST_AUTO_TEST_CASE(EdgeHandle)
{
    const auto& edges = handle.edges();
    BOOST_TEST(std::distance(edges.begin(), edges.end()) == edge_cnt);
    BOOST_TEST(edges.size() == edge_cnt);
    std::set<int> vs;
    for(const auto & e : handle.edges()) {
        vs.insert(e.vertexID(0));
        vs.insert(e.vertexID(1));
    }
    BOOST_TEST(vs.size() == vertex_cnt);
}

BOOST_AUTO_TEST_CASE(TriangleHandle)
{
    const auto& triangles = handle.triangles();
    BOOST_TEST(std::distance(triangles.begin(), triangles.end()) == triangle_cnt);
    BOOST_TEST(triangles.size() == triangle_cnt);
    std::set<int> vs;
    for(const auto & t : handle.triangles()) {
        vs.insert(t.vertexID(0));
        vs.insert(t.vertexID(1));
        vs.insert(t.vertexID(2));
    }
    BOOST_TEST(vs.size() == vertex_cnt);
}

BOOST_AUTO_TEST_SUITE_END()
