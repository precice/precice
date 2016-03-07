// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "SpacetreeTestScenarios.hpp"
#include "spacetree/ExportSpacetree.hpp"
#include "io/ExportVTK.hpp"
#include "utils/Helpers.hpp"
#include "query/FindClosest.hpp"
#include "query/FindVoxelContent.hpp"
#include "geometry/Sphere.hpp"
#include "geometry/Cuboid.hpp"
#include "mesh/Merge.hpp"

namespace precice {
namespace spacetree {
namespace tests {

tarch::logging::Log SpacetreeTestScenarios::
  _log("precice::spacetree::tests::SpacetreeTestScenarios");

SpacetreeTestScenarios:: SpacetreeTestScenarios
(
  const std::string& testName,
  SpacetreeFactory&  factory )
:
  _testName(testName),
  _factory(factory)
{}

void SpacetreeTestScenarios:: testSearchPosition()
{
  preciceTrace ( "testSearchPosition()" );
  typedef Spacetree Tree;
  int min = Tree::minElementsToRefineCell;
  Tree::minElementsToRefineCell = 1;
  int dim = 2;
  { // 2D
    using utils::Vector2D;
    // Create mesh
    bool flipNormals = false;
    mesh::PtrMesh mesh(new mesh::Mesh("TestMesh", dim, flipNormals));
    mesh::Vertex& v0 = mesh->createVertex(Vector2D(0.0, 0.0));
    mesh::Vertex& v00 = mesh->createVertex(Vector2D(0.3, 0.0));
    mesh::Vertex& v01 = mesh->createVertex(Vector2D(0.7, 0.0));
    mesh::Vertex& v1 = mesh->createVertex(Vector2D(1.0, 0.0));
    mesh::Vertex& v2 = mesh->createVertex(Vector2D(1.0, 1.0));
    mesh::Vertex& v3 = mesh->createVertex(Vector2D(0.0, 1.0));
    mesh::Vertex& v30 = mesh->createVertex(Vector2D(0.0, 0.85));
    mesh::Vertex& v31 = mesh->createVertex(Vector2D(0.0, 0.7));
    mesh::Vertex& v32 = mesh->createVertex(Vector2D(0.0, 0.3));
    mesh->createEdge(v0, v00);
    mesh->createEdge(v00, v01);
    mesh->createEdge(v01, v1);
    mesh->createEdge(v1, v2);
    mesh->createEdge(v2, v3);
    mesh->createEdge(v3, v30);
    mesh->createEdge(v30, v31);
    mesh->createEdge(v31, v32);
    mesh->createEdge(v32, v0);
    mesh->computeState();

    // Create and initialize spacetree
    utils::DynVector center(Vector2D(0.5, 0.5));
    utils::DynVector halflengths(2, 2.0);
    double upperRefinementLimit = 0.06125;
    PtrSpacetree spacetree = _factory.createSpacetree(center, halflengths,
                                                      upperRefinementLimit);
    spacetree->addMesh(mesh);
    spacetree->initialize();

#   ifdef Debug
    int testNumber = 1;
#   endif

    // Perform tests
    // Outside positions
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 1.5, -0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 1.5, 0.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 1.5, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 1.5, 1.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 1.5, 1.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( -0.5, -0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( -0.5, 0.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( -0.5, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( -0.5, 1.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( -0.5, 1.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 0.5, 1.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 0.5, -0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }

    // Outside eps positions
    double eps = tarch::la::NUMERICAL_ZERO_DIFFERENCE;
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 1.0 + 10.0 * eps, -0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 1.0 + 10.0 * eps, 0.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 1.0 + 10.0 * eps, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 1.0 + 10.0 * eps, 1.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 1.0 + 10.0 * eps, 1.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( -10.0 * eps, -0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( -10.0 * eps, 0.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( -10.0 * eps, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( -10.0 * eps, 1.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( -10.0 * eps, 1.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 0.5, 1.0 + 10.0 * eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 0.5, -0.5 - 10.0 * eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }

    // Touching
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 0.0, 0.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 0.5, 0.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 1.0, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 1.0, 1.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 0.5, 1.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 0.0, 1.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 0.0, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }

    // Touching eps
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 0.0 + 0.1 * eps, 0.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 0.5, 0.0 - 0.1 * eps);
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 1.0 + 0.1 * eps, 0.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 1.0 - 0.1 * eps, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 1.0 + 0.1 * eps, 1.0 + 0.1 * eps);
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 0.5, 1.0 - 0.1 * eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 0.0 - 0.1 * eps, 1.0 + 0.1 * eps);
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 0.0 - 0.1 * eps, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }

    // Inside eps
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 0.0 + 10.0 * eps, 0.0 + 10.0 * eps);
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 0.5, 0.0 + 10.0 * eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 1.0 - 10.0 * eps, 0.0 + 10.0 * eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 1.0 - 10.0 * eps, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 1.0 - 10.0 * eps, 1.0 - 10.0 * eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 0.5, 1.0 - 10.0 * eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 0.0 + 10.0 * eps, 1.0 - 10.0 * eps);
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 0.0 + 10.0 * eps, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }

    // Inside
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 0.5, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 0.1, 0.1 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 0.9, 0.9 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 0.1, 0.9 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  2d test "  << testNumber++);
      Vector2D point ( 0.9, 0.1 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }

    // Export geometric results to vtk files
    ExportSpacetree exportSpacetree ( _testName + "-testSearchPosition-tree-2d" );
    exportSpacetree.doExport(*spacetree);
    io::ExportVTK exportMesh(true);
    exportMesh.doExport ( _testName + "-testSearchPosition-mesh-2d", *mesh );
  }

  dim = 3;
  {
    using utils::Vector3D;
    // Create mesh
    bool flipNormals = false;
    mesh::PtrMesh mesh(new mesh::Mesh("TestMesh", dim, flipNormals));
    mesh::Vertex& v000 = mesh->createVertex ( Vector3D(0.0, 0.0, 0.0) );
    mesh::Vertex& v001 = mesh->createVertex ( Vector3D(0.0, 0.0, 1.0) );
    mesh::Vertex& v010 = mesh->createVertex ( Vector3D(0.0, 1.0, 0.0) );
    mesh::Vertex& v011 = mesh->createVertex ( Vector3D(0.0, 1.0, 1.0) );
    mesh::Vertex& v100 = mesh->createVertex ( Vector3D(1.0, 0.0, 0.0) );
    mesh::Vertex& v101 = mesh->createVertex ( Vector3D(1.0, 0.0, 1.0) );
    mesh::Vertex& v110 = mesh->createVertex ( Vector3D(1.0, 1.0, 0.0) );
    mesh::Vertex& v111 = mesh->createVertex ( Vector3D(1.0, 1.0, 1.0) );

    mesh::Edge& e000to100 = mesh->createEdge ( v000, v100 );
    mesh::Edge& e010to110 = mesh->createEdge ( v010, v110 );
    mesh::Edge& e001to101 = mesh->createEdge ( v001, v101 );
    mesh::Edge& e011to111 = mesh->createEdge ( v011, v111 );

    mesh::Edge& e000to010 = mesh->createEdge ( v000, v010 );
    mesh::Edge& e100to110 = mesh->createEdge ( v100, v110 );
    mesh::Edge& e001to011 = mesh->createEdge ( v001, v011 );
    mesh::Edge& e101to111 = mesh->createEdge ( v101, v111 );

    mesh::Edge& e000to001 = mesh->createEdge ( v000, v001 );
    mesh::Edge& e100to101 = mesh->createEdge ( v100, v101 );
    mesh::Edge& e010to011 = mesh->createEdge ( v010, v011 );
    mesh::Edge& e110to111 = mesh->createEdge ( v110, v111 );

    mesh::Edge& e000to011 = mesh->createEdge ( v000, v011 );
    mesh::Edge& e100to111 = mesh->createEdge ( v100, v111 );
    mesh::Edge& e100to001 = mesh->createEdge ( v100, v001 );
    mesh::Edge& e110to011 = mesh->createEdge ( v110, v011 );
    mesh::Edge& e000to110 = mesh->createEdge ( v000, v110 );
    mesh::Edge& e001to111 = mesh->createEdge ( v001, v111 );

    mesh->createTriangle ( e000to001, e001to011, e000to011 ); // x = 0
    mesh->createTriangle ( e000to010, e000to011, e010to011 );
    mesh->createTriangle ( e100to101, e100to111, e101to111 ); // x = 1
    mesh->createTriangle ( e100to110, e110to111, e100to111 );

    mesh->createTriangle ( e000to100, e100to001, e000to001 ); // y = 0
    mesh->createTriangle ( e100to101, e001to101, e100to001 );
    mesh->createTriangle ( e010to110, e010to011, e110to011 ); // y = 1
    mesh->createTriangle ( e110to111, e110to011, e011to111 );

    mesh->createTriangle ( e000to100, e000to110, e100to110 ); // z = 0
    mesh->createTriangle ( e000to010, e010to110, e000to110 );
    mesh->createTriangle ( e001to101, e101to111, e001to111 ); // z = 1
    mesh->createTriangle ( e001to011, e001to111, e011to111 );

    mesh->computeState();

    // Create and initialize spacetree
    utils::DynVector center(dim, 0.5);
    utils::DynVector halflengths(dim, 2.0);
    double upperRefinementLimit = 0.06125;
    PtrSpacetree spacetree = _factory.createSpacetree( center, halflengths,
                                                       upperRefinementLimit );
    spacetree->addMesh(mesh);
    spacetree->initialize();

#   ifdef Debug
    int testNumber = 1;
#   endif

    // Perform tests
    // Outside
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.5, 1.5, 1.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.5, 1.5, -0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.5, -0.5, 1.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.5, -0.5, -0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( -0.5, 1.5, 1.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( -0.5, 1.5, -0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( -0.5, -0.5, 1.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( -0.5, -0.5, -0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.5, 0.5, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( -0.5, 0.5, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.5, 1.5, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.5, -0.5, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.5, 0.5, 1.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.5, 0.5, -0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }

    // Outside eps
    double eps = tarch::la::NUMERICAL_ZERO_DIFFERENCE;
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.0 + 10.0*eps, 1.0 + 10.0*eps, 1.0 + 10.0*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.0 + 10.0*eps, 1.0 + 10.0*eps, -10.0*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.0 + 10.0*eps, -10.0*eps, 1.0 + 10.0*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.0 + 10.0*eps, -10.0*eps, -10.0*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( -10.0*eps, 1.0 + 10.0*eps, 1.0 + 10.0*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( -10.0*eps, 1.0 + 10.0*eps, -10.0*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( -10.0*eps, -10.0*eps, 1.0 + 10.0*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( -10.0*eps, -10.0*eps, -10.0*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.0 + 10.0*eps, 0.5, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( -10.0*eps, 0.5, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.5, 1.0 + 10.0*eps, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.5, -10.0*eps, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.5, 0.5, 1.0 + 10.0*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.5, 0.5, -10.0*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOutsideOfGeometry() );
    }

    // Touching eps
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.0 + 0.1*eps, 1.0 + 0.1*eps, 1.0 + 0.1*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.0 + 0.1*eps, 1.0 + 0.1*eps, -0.1*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.0 + 0.1*eps, -0.1*eps, 1.0 + 0.1*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.0 + 0.1*eps, -0.1*eps, -0.1*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( -0.1*eps, 1.0 + 0.1*eps, 1.0 + 0.1*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( -0.1*eps, 1.0 + 0.1*eps, -0.1*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( -0.1*eps, -0.1*eps, 1.0 + 0.1*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( -0.1*eps, -0.1*eps, -0.1*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.0 + 0.1*eps, 0.5, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( -0.1*eps, 0.5, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.5, 1.0 + 0.1*eps, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.5, -0.1*eps, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.5, 0.5, 1.0 + 0.1*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.5, 0.5, -0.1*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }

    // Touching
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.0, 1.0, 1.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.0, 1.0, 0.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.0, 0.0, 1.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.0, 0.0, 0.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.0, 1.0, 1.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.0, 1.0, 0.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.0, 0.0, 1.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.0, 0.0, 0.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.0, 0.5, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.0, 0.5, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.5, 1.0, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.5, 0.0, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.5, 0.5, 1.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.5, 0.5, 0.0 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionOnGeometry() );
    }

    // Inside eps
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.0 - 10.0*eps, 1.0 - 10.0*eps, 1.0 - 10.0*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.0 - 10.0*eps, 1.0 - 10.0*eps, 10.0*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.0 - 10.0*eps, 10.0*eps, 1.0 - 10.0*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.0 - 10.0*eps, 10.0*eps, 10.0*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 10.0*eps, 1.0 - 10.0*eps, 1.0 - 10.0*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 10.0*eps, 1.0 - 10.0*eps, 10.0*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 10.0*eps, 10.0*eps, 1.0 - 10.0*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 10.0*eps, 10.0*eps, 10.0*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 1.0 - 10.0*eps, 0.5, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 10.0*eps, 0.5, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.5, 1.0 - 10.0*eps, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.5, 10.0*eps, 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.5, 0.5, 1.0 - 10.0*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.5, 0.5, 10.0*eps );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }

    // Inside
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.5 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.1 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.9 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.1, 0.9, 0.9 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.1, 0.1, 0.9 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.9, 0.1, 0.1 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.9, 0.9, 0.1 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }
    {
      preciceDebug("  3d test "  << testNumber++);
      Vector3D point ( 0.9, 0.1, 0.9 );
      int pos = spacetree->searchPosition ( point );
      validateEquals ( pos, Tree::positionInsideOfGeometry() );
    }

    // Export geometric results to vtk files
    ExportSpacetree exportSpacetree ( _testName + "-testSearchPosition-tree-3d" );
    exportSpacetree.doExport(*spacetree);
    io::ExportVTK exportMesh(true);
    exportMesh.doExport ( _testName + "-testSearchPosition-mesh-3d", *mesh );
  }
  Tree::minElementsToRefineCell = min;
}

void SpacetreeTestScenarios:: testSearchDistance()
{
  preciceTrace ( "testSearchDistance()" );
  int min = Spacetree::minElementsToRefineCell;
  Spacetree::minElementsToRefineCell = 1;
  int dim = 2;
  { // 2D
    using utils::Vector2D;
    // Create mesh
    bool flipNormals = false;
    mesh::PtrMesh mesh(new mesh::Mesh("TestMesh", dim, flipNormals));
    mesh::Vertex& v0 = mesh->createVertex ( Vector2D(0.0, 0.0) );
    mesh::Vertex& v00 = mesh->createVertex ( Vector2D(0.3, 0.0) );
    mesh::Vertex& v01 = mesh->createVertex ( Vector2D(0.7, 0.0) );
    mesh::Vertex& v1 = mesh->createVertex ( Vector2D(1.0, 0.0) );
    mesh::Vertex& v2 = mesh->createVertex ( Vector2D(1.0, 1.0) );
    mesh::Vertex& v3 = mesh->createVertex ( Vector2D(0.0, 1.0) );
    mesh::Vertex& v30 = mesh->createVertex ( Vector2D(0.0, 0.85) );
    mesh::Vertex& v31 = mesh->createVertex ( Vector2D(0.0, 0.7) );
    mesh::Vertex& v32 = mesh->createVertex ( Vector2D(0.0, 0.3) );
    mesh->createEdge ( v0, v00 );
    mesh->createEdge ( v00, v01 );
    mesh->createEdge ( v01, v1 );
    mesh->createEdge ( v1, v2 );
    mesh->createEdge ( v2, v3 );
    mesh->createEdge ( v3, v30 );
    mesh->createEdge ( v30, v31 );
    mesh->createEdge ( v31, v32 );
    mesh->createEdge ( v32, v0 );
    mesh->computeState();

    // Create and initialize spacetree
    utils::DynVector center(Vector2D(0.5, 0.5));
    utils::DynVector halflengths(2, 1.0);
    double upperRefinementLimit = 0.06125;
    PtrSpacetree spacetree = _factory.createSpacetree ( center, halflengths,
                                                        upperRefinementLimit );
    spacetree->addMesh(mesh);
    spacetree->initialize();

    // Perform tests
    //query::ExportVTKNeighbors exportNeighbors;
    {
      Vector2D pos ( 0.3, 0.3 );
      query::FindClosest search ( pos );
      spacetree->searchDistance ( search );
      query::ClosestElement closest = search.getClosest ();
      validateNumericalEquals ( closest.distance, -0.3 );
      //exportNeighbors.addNeighbors ( pos, closest );
    }
    {
      Vector2D pos ( 0.5, 0.5 );
      query::FindClosest search ( pos );
      spacetree->searchDistance ( search );
      query::ClosestElement closest = search.getClosest ();
      validateNumericalEquals ( closest.distance, -0.5 );
      //exportNeighbors.addNeighbors ( pos, closest );
    }
    {
      Vector2D pos ( 0.7, 0.7 );
      query::FindClosest search ( pos );
      spacetree->searchDistance ( search );
      query::ClosestElement closest = search.getClosest ();
      validateNumericalEquals ( closest.distance, -0.3 );
      //exportNeighbors.addNeighbors ( pos, closest );
    }
    {
      Vector2D pos ( -0.1, 1.0 );
      query::FindClosest search ( pos );
      spacetree->searchDistance ( search );
      query::ClosestElement closest = search.getClosest ();
      validateNumericalEquals ( closest.distance, 0.1 );
      //exportNeighbors.addNeighbors ( pos, closest );
    }
    {
      Vector2D pos ( 1.0, 1.1 );
      query::FindClosest search ( pos );
      spacetree->searchDistance ( search );
      query::ClosestElement closest = search.getClosest ();
      validateNumericalEquals ( closest.distance, 0.1 );
      //exportNeighbors.addNeighbors ( pos, closest );
    }
    {
      Vector2D pos ( 0.5, 1.2 );
      query::FindClosest search ( pos );
      spacetree->searchDistance ( search );
      query::ClosestElement closest = search.getClosest ();
      validateNumericalEquals ( closest.distance, 0.2 );
      //exportNeighbors.addNeighbors ( pos, closest );
    }
    {
      Vector2D pos ( -0.4, 0.7 );
      query::FindClosest search ( pos );
      spacetree->searchDistance ( search );
      query::ClosestElement closest = search.getClosest ();
      validateNumericalEquals ( closest.distance, 0.4 );
      //exportNeighbors.addNeighbors ( pos, closest );
    }
    {
      Vector2D pos ( 0.7, -0.5 );
      query::FindClosest search ( pos );
      spacetree->searchDistance ( search );
      query::ClosestElement closest = search.getClosest ();
      validateNumericalEquals ( closest.distance, 0.5 );
      //exportNeighbors.addNeighbors ( pos, closest );
    }

    // Export geometric results to vtk files
    ExportSpacetree exportSpacetree ( _testName + "-testSearchDistance-tree-2D" );
    exportSpacetree.doExport(*spacetree);
    io::ExportVTK exportMesh(true);
    exportMesh.doExport ( _testName + "-testSearchDistance-mesh-2D", *mesh );
//    exportNeighbors.exportNeighbors ( "SpacetreeTestScenarios-testSearchDistance-neighbors-2D.vtk");
  }

  dim = 3;
  { // 3D
    using utils::Vector3D;
    // Create mesh
    bool flipNormals = false;
    mesh::PtrMesh mesh(new mesh::Mesh("TestMesh", dim, flipNormals));
    mesh::Vertex& v000 = mesh->createVertex ( Vector3D(0.0, 0.0, 0.0) );
    mesh::Vertex& v001 = mesh->createVertex ( Vector3D(0.0, 0.0, 1.0) );
    mesh::Vertex& v010 = mesh->createVertex ( Vector3D(0.0, 1.0, 0.0) );
    mesh::Vertex& v011 = mesh->createVertex ( Vector3D(0.0, 1.0, 1.0) );
    mesh::Vertex& v100 = mesh->createVertex ( Vector3D(1.0, 0.0, 0.0) );
    mesh::Vertex& v101 = mesh->createVertex ( Vector3D(1.0, 0.0, 1.0) );
    mesh::Vertex& v110 = mesh->createVertex ( Vector3D(1.0, 1.0, 0.0) );
    mesh::Vertex& v111 = mesh->createVertex ( Vector3D(1.0, 1.0, 1.0) );

    mesh::Edge& e000to100 = mesh->createEdge ( v000, v100 );
    mesh::Edge& e010to110 = mesh->createEdge ( v010, v110 );
    mesh::Edge& e001to101 = mesh->createEdge ( v001, v101 );
    mesh::Edge& e011to111 = mesh->createEdge ( v011, v111 );

    mesh::Edge& e000to010 = mesh->createEdge ( v000, v010 );
    mesh::Edge& e100to110 = mesh->createEdge ( v100, v110 );
    mesh::Edge& e001to011 = mesh->createEdge ( v001, v011 );
    mesh::Edge& e101to111 = mesh->createEdge ( v101, v111 );

    mesh::Edge& e000to001 = mesh->createEdge ( v000, v001 );
    mesh::Edge& e100to101 = mesh->createEdge ( v100, v101 );
    mesh::Edge& e010to011 = mesh->createEdge ( v010, v011 );
    mesh::Edge& e110to111 = mesh->createEdge ( v110, v111 );

    mesh::Edge& e000to011 = mesh->createEdge ( v000, v011 );
    mesh::Edge& e100to111 = mesh->createEdge ( v100, v111 );
    mesh::Edge& e100to001 = mesh->createEdge ( v100, v001 );
    mesh::Edge& e110to011 = mesh->createEdge ( v110, v011 );
    mesh::Edge& e000to110 = mesh->createEdge ( v000, v110 );
    mesh::Edge& e001to111 = mesh->createEdge ( v001, v111 );

    mesh->createTriangle ( e000to001, e001to011, e000to011 ); // x = 0
    mesh->createTriangle ( e000to010, e000to011, e010to011 );
    mesh->createTriangle ( e100to101, e100to111, e101to111 ); // x = 1
    mesh->createTriangle ( e100to110, e110to111, e100to111 );

    mesh->createTriangle ( e000to100, e100to001, e000to001 ); // y = 0
    mesh->createTriangle ( e100to101, e001to101, e100to001 );
    mesh->createTriangle ( e010to110, e010to011, e110to011 ); // y = 1
    mesh->createTriangle ( e110to111, e110to011, e011to111 );

    mesh->createTriangle ( e000to100, e000to110, e100to110 ); // z = 0
    mesh->createTriangle ( e000to010, e010to110, e000to110 );
    mesh->createTriangle ( e001to101, e101to111, e001to111 ); // z = 1
    mesh->createTriangle ( e001to011, e001to111, e011to111 );

    mesh->computeState();

    // Create and initialize spacetree
    utils::DynVector center(dim, 0.5);
    utils::DynVector halflengths(dim, 1.0);
    double upperRefinementLimit = 0.06125;
    PtrSpacetree spacetree = _factory.createSpacetree ( center, halflengths,
                                                        upperRefinementLimit );
    spacetree->addMesh(mesh);
    spacetree->initialize();

    // Perform tests
//    query::ExportVTKNeighbors exportNeighbors;
    {
      Vector3D pos ( 0.3 );
      query::FindClosest search ( pos );
      spacetree->searchDistance ( search );
      query::ClosestElement closest = search.getClosest ();
      validateNumericalEquals ( closest.distance, -0.3 );
      //exportNeighbors.addNeighbors ( pos, closest );
    }
    {
      Vector3D pos ( 0.7 );
      query::FindClosest search ( pos );
      spacetree->searchDistance ( search );
      query::ClosestElement closest = search.getClosest ();
      validateNumericalEquals ( closest.distance, -0.3 );
      //exportNeighbors.addNeighbors ( pos, closest );
    }
    {
      Vector3D pos ( 0.5 );
      query::FindClosest search ( pos );
      spacetree->searchDistance ( search );
      query::ClosestElement closest = search.getClosest ();
      validateNumericalEquals ( closest.distance, -0.5 );
      //exportNeighbors.addNeighbors ( pos, closest );
    }
    {
      Vector3D pos ( 0.5, 0.5, 0.8 );
      query::FindClosest search ( pos );
      spacetree->searchDistance ( search );
      query::ClosestElement closest = search.getClosest ();
      validateNumericalEquals ( closest.distance, -0.2 );
      //exportNeighbors.addNeighbors ( pos, closest );
    }
    {
      Vector3D pos ( 0.5, 0.8, 0.5 );
      query::FindClosest search ( pos );
      spacetree->searchDistance ( search );
      query::ClosestElement closest = search.getClosest ();
      validateNumericalEquals ( closest.distance, -0.2 );
      //exportNeighbors.addNeighbors ( pos, closest );
    }
    {
      Vector3D pos ( 0.8, 0.5, 0.5 );
      query::FindClosest search ( pos );
      spacetree->searchDistance ( search );
      query::ClosestElement closest = search.getClosest ();
      validateNumericalEquals ( closest.distance, -0.2 );
      //exportNeighbors.addNeighbors ( pos, closest );
    }

    // Export geometric results to vtk files
    ExportSpacetree exportSpacetree ( _testName + "-testSearchDistance-tree-3D" );
    exportSpacetree.doExport ( *spacetree );
    io::ExportVTK exportMesh(true);
    exportMesh.doExport ( _testName + "-testSearchDistance-mesh-3D", *mesh );
//    exportNeighbors.exportNeighbors ( "SpacetreeTestScenarios-testSearchDistance-neighbors-3D.vtk");
  }
  Spacetree::minElementsToRefineCell = min;
}

void SpacetreeTestScenarios:: testNeighborSearch()
{
  preciceTrace ( "testNeighborSearch()" );
  using utils::Vector2D;
  int dim = 2;
  double sphereRadius = 3.0;
  bool flipNormals = true;
  mesh::PtrMesh mesh(new mesh::Mesh("test-sphere", dim, flipNormals));
  geometry::Sphere(utils::DynVector(dim,0.0), 0.01, sphereRadius).create(*mesh);

  io::ExportVTK exportVTK(true);
  exportVTK.doExport ( _testName + "-testNeighborSearch-sphere.vtk", *mesh );

  utils::DynVector center(dim, 0.0);
  utils::DynVector h(dim, 10.0);
  PtrSpacetree spacetree = _factory.createSpacetree ( center, h, 0.05 );
  spacetree->addMesh(mesh);
  spacetree->initialize();
  //query::ExportVTKNeighbors exportNeighborsSpacetree;
  //query::ExportVTKNeighbors exportNeighbors;

  Vector2D point ( 0.0 );
  Vector2D point2 ( sphereRadius + sphereRadius / 20.0, 0.0 );
  Vector2D point3 ( sphereRadius, sphereRadius );
  Vector2D point4 ( sphereRadius * 0.8, sphereRadius * -0.85 );
  Vector2D point5 ( sphereRadius * 0.1, sphereRadius * -0.5 );
  double width = 6.0 * sphereRadius;
  int steps = 20;
  double stepWidth = width/(double)steps;
  for ( int i=0; i < steps; i++ ) {
    Vector2D point_i (-width/2.0 + (double)i*stepWidth, sphereRadius * 1.5);
    query::FindClosest findNeighbors0 ( point_i );
    findNeighbors0 ( *mesh );
    //exportNeighbors.addNeighbors ( point_i, findNeighbors0.getClosest() );
    findNeighbors0.reset ();
    spacetree->searchDistance ( findNeighbors0 );
    //exportNeighborsspacetree->addNeighbors ( point_i, findNeighbors0.getClosest() );
    assignList(point_i) = -width/2.0 + (double)i*stepWidth, - sphereRadius * 0.8;
    query::FindClosest findNeighbors1 ( point_i );
    findNeighbors1 ( *mesh );
    //exportNeighbors.addNeighbors ( point_i, findNeighbors1.getClosest() );
    findNeighbors1.reset ();
    spacetree->searchDistance ( findNeighbors1 );
    //exportNeighborsspacetree->addNeighbors ( point_i, findNeighbors1.getClosest() );
    assignList(point_i) = -width/2.0 + (double)i*stepWidth, - sphereRadius * 0.2;
    query::FindClosest findNeighbors2 ( point_i );
    findNeighbors2 ( *mesh );
    //exportNeighbors.addNeighbors ( point_i, findNeighbors2.getClosest() );
    findNeighbors2.reset ();
    spacetree->searchDistance ( findNeighbors2 );
    //exportNeighborsspacetree->addNeighbors ( point_i, findNeighbors2.getClosest() );
  }

  query::FindClosest findNeighbors0 ( point );
  findNeighbors0 ( *mesh );
  //exportNeighbors.addNeighbors ( point, findNeighbors0.getClosest() );
  findNeighbors0.reset();
  spacetree->searchDistance ( findNeighbors0 );
  //exportNeighborsspacetree->addNeighbors ( point, findNeighbors0.getClosest() );

  query::FindClosest findNeighbors1 ( point2 );
  findNeighbors1 ( *mesh );
  //exportNeighbors.addNeighbors ( point2, findNeighbors1.getClosest() );
  findNeighbors1.reset ();
  spacetree->searchDistance ( findNeighbors1 );
  //exportNeighborsspacetree->addNeighbors ( point2, findNeighbors1.getClosest() );

  query::FindClosest findNeighbors2 ( point3 );
  findNeighbors2 ( *mesh );
  //exportNeighbors.addNeighbors ( point3, findNeighbors2.getClosest() );
  findNeighbors2.reset ();
  spacetree->searchDistance ( findNeighbors2 );
  //exportNeighborsspacetree->addNeighbors ( point3, findNeighbors2.getClosest() );

  query::FindClosest findNeighbors3 ( point4 );
  findNeighbors3 ( *mesh );
  //exportNeighbors.addNeighbors ( point4, findNeighbors3.getClosest() );
  findNeighbors3.reset ();
  spacetree->searchDistance ( findNeighbors3 );
  //exportNeighborsspacetree->addNeighbors ( point4, findNeighbors3.getClosest() );

  query::FindClosest findNeighbors4 ( point5 );
  findNeighbors4 ( *mesh );
  //exportNeighbors.addNeighbors ( point5, findNeighbors4.getClosest() );
  findNeighbors4.reset ();
  spacetree->searchDistance ( findNeighbors4 );
  //exportNeighborsspacetree->addNeighbors ( point5, findNeighbors4.getClosest() );

  ExportSpacetree exportSpacetree (_testName + "-testNeighborSearch-spacetree.vtk");
  exportSpacetree.doExport ( *spacetree );
//  exportNeighbors.exportNeighbors (
//      "SpacetreeTestScenarios_testNeighborSearch_neighbors.vtk");
//  exportNeighborsspacetree->exportNeighbors (
//      "SpacetreeTestScenarios_testNeighborSearch_neighbors_withspacetree->vtk");
}

void SpacetreeTestScenarios:: testSearchContentVertices()
{
  preciceTrace ( "testSearchContentVertices()" );
  for ( int dim=2; dim <= 3; dim++ ){
    for ( int testDim=0; testDim < dim; testDim++ ) {
      bool positiveDirection = true;
      bool negativeDirection = false;
      utils::DynVector offset ( dim, 0.0 );
      performTestSearchContentVertices ( testDim, positiveDirection, offset );
      performTestSearchContentVertices ( testDim, negativeDirection, offset );
      assign(offset) = -1.0 + 10.0 * tarch::la::NUMERICAL_ZERO_DIFFERENCE;
      performTestSearchContentVertices ( testDim, positiveDirection, offset );
      performTestSearchContentVertices ( testDim, negativeDirection, offset );
      assign(offset) = 1.0 - 10.0 * tarch::la::NUMERICAL_ZERO_DIFFERENCE;
      performTestSearchContentVertices ( testDim, positiveDirection, offset );
      performTestSearchContentVertices ( testDim, negativeDirection, offset );
    }
  }
}

void SpacetreeTestScenarios:: performTestSearchContentVertices
(
  int                     testDim,
  bool                    positive,
  const utils::DynVector& offset )
{
  preciceTrace3 ( "performTestSearchContentVertices()", testDim, positive, offset );
  int min = Spacetree::minElementsToRefineCell;
  Spacetree::minElementsToRefineCell = 1;
  int dim = offset.size();
  assertion ( not tarch::la::oneGreater(offset, utils::DynVector(dim,1.0)) );
  assertion ( tarch::la::allGreater(offset, utils::DynVector(dim,-1.0)) );
  bool flipNormals = false;
  mesh::PtrMesh mesh(new mesh::Mesh("TestMesh", dim, flipNormals));
  utils::DynVector coords(offset);
  mesh::Vertex& vertex = mesh->createVertex(coords);

  utils::DynVector center(dim, 0.0);
  utils::DynVector halflengths(dim, 1.0);
  query::FindVoxelContent::BoundaryInclusion includeBounds =
      query::FindVoxelContent::INCLUDE_BOUNDARY;
  query::FindVoxelContent::BoundaryInclusion excludeBounds =
      query::FindVoxelContent::EXCLUDE_BOUNDARY;
  query::FindVoxelContent findIncluded ( center, halflengths, includeBounds );
  query::FindVoxelContent findExcluded ( center, halflengths, excludeBounds );

  utils::DynVector treeOffset(dim, 0.0);
  utils::DynVector treeHalflengths(dim, 2.0);
  std::vector<double> refinementLimits;
  refinementLimits += 2.0, 1.0, 0.5, 0.25;
  std::vector<PtrSpacetree> treesInc;
  std::vector<PtrSpacetree> treesExc;
  for ( double limit : refinementLimits ) {
    treesInc +=_factory.createSpacetree(treeOffset, treeHalflengths, limit);
    treesInc.back()->addMesh(mesh);
    treesExc += _factory.createSpacetree(treeOffset, treeHalflengths, limit);
    treesExc.back()->addMesh(mesh);
  }

  assertion ( testDim >= 0 );
  assertion ( testDim < dim );

  double sign = positive ? 1.0 : -1.0;
  int size = 0;
  mesh::Merge merge;

  // Outside
  coords[testDim] = sign * 2.0;
  vertex.setCoords ( coords );
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ){
    preciceDebug ( "outside i= " << i );
    treesInc[i]->searchContent(findIncluded);
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().vertices().size();
    validateEquals ( size, 0 );
    size = findExcluded.content().vertices().size();
    validateEquals ( size, 0 );
  }

  // Outside eps
  coords[testDim] = sign * (1.0 + 10.0 * tarch::la::NUMERICAL_ZERO_DIFFERENCE);
  vertex.setCoords(coords);
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ){
    preciceDebug ( "outside i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().vertices().size();
    validateEquals ( size, 0 );
    size = findExcluded.content().vertices().size();
    validateEquals ( size, 0 );
  }

  // Outside eps
  coords[testDim] = sign * (1.0 + 10.0 * tarch::la::NUMERICAL_ZERO_DIFFERENCE);
  vertex.setCoords(coords);
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ){
    preciceDebug ( "outside i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().vertices().size();
    validateEquals ( size, 0 );
    size = findExcluded.content().vertices().size();
    validateEquals ( size, 0 );
  }

  // Touching + eps
  coords[testDim] = sign * (1.0 + tarch::la::NUMERICAL_ZERO_DIFFERENCE);
  vertex.setCoords(coords);
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ){
    preciceDebug ( "outside i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().vertices().size();
    validateEquals ( size, 1 );
    findIncluded.clear ();
    size = findExcluded.content().vertices().size();
    validateEquals ( size, 0 );
  }

  // Touching
  coords[testDim] = sign * 1.0;
  vertex.setCoords(coords);
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ){
    preciceDebug ( "outside i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().vertices().size();
    validateEquals ( size, 1 );
    findIncluded.clear ();
    size = findExcluded.content().vertices().size();
    validateEquals ( size, 0 );
  }

  // Touching - eps
  coords[testDim] = sign * (1.0 - tarch::la::NUMERICAL_ZERO_DIFFERENCE);
  vertex.setCoords(coords);
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ){
    preciceDebug ( "outside i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().vertices().size();
    validateEquals ( size, 1 );
    findIncluded.clear ();
    size = findExcluded.content().vertices().size();
    validateEquals ( size, 0 );
  }

  // Inside eps
  coords[testDim] = sign * (1.0 - 10.0 * tarch::la::NUMERICAL_ZERO_DIFFERENCE);
  vertex.setCoords(coords);
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ){
    preciceDebug ( "outside i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().vertices().size();
    validateEquals ( size, 1 );
    findIncluded.clear ();
    size = findExcluded.content().vertices().size();
    validateEquals ( size, 1 );
    findExcluded.clear ();
  }

  // Inside
  coords[testDim] = sign * 0.9;
  vertex.setCoords(coords);
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ){
    preciceDebug ( "outside i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().vertices().size();
    validateEquals ( size, 1 );
    findIncluded.clear ();
    size = findExcluded.content().vertices().size();
    validateEquals ( size, 1 );
    findExcluded.clear ();
  }
  Spacetree::minElementsToRefineCell = min;
}

void SpacetreeTestScenarios:: testSearchContentEdges()
{
  preciceTrace("testSearchContentEdges()");
  for (int dim=2; dim <= 3; dim++){
    for (int testDim=0; testDim < dim; testDim++){
      bool positiveDirection = true;
      bool negativeDirection = false;
      utils::DynVector offset(dim, 0.0);
      performTestSearchContentEdges(testDim, positiveDirection, offset);
      performTestSearchContentEdges(testDim, negativeDirection, offset);
      assign(offset) = -1.0 + 10.0 * tarch::la::NUMERICAL_ZERO_DIFFERENCE;
      performTestSearchContentEdges(testDim, positiveDirection, offset);
      performTestSearchContentEdges(testDim, negativeDirection, offset);
      assign(offset) = 1.0 - 10.0 * tarch::la::NUMERICAL_ZERO_DIFFERENCE;
      performTestSearchContentEdges(testDim, positiveDirection, offset);
      performTestSearchContentEdges(testDim, negativeDirection, offset);
    }
  }

  // Special test
  preciceDebug("Reproduce bug test");
  Spacetree::minElementsToRefineCell = 1;
  int dim = 3;
  bool flipNormals = false;
  mesh::PtrMesh mesh(new mesh::Mesh("TestMesh", dim, flipNormals));
  utils::DynVector coords0(dim), coords1(dim), coords2(dim);
  coords0 = 1.0, 0.0, 1.0;
  coords1 = 0.0, 1.0, 1.0;
  coords2 = 1.0, 1.0, 1.0;
  mesh::Vertex& v0 = mesh->createVertex(coords0);
  mesh::Vertex& v1 = mesh->createVertex(coords1);
  mesh::Vertex& v2 = mesh->createVertex(coords2);
  mesh::Edge& e0 = mesh->createEdge(v0, v1);
  mesh::Edge& e1 = mesh->createEdge(v1, v2);
  mesh::Edge& e2 = mesh->createEdge(v0, v2);
  mesh->createTriangle(e2, e1, e0);
  mesh->computeState();

  utils::DynVector center(dim), halflengths(dim);
  center = 0.5, 0.5, 0.5;
  halflengths = center;
  PtrSpacetree tree = _factory.createSpacetree(center, halflengths, 1.0/8.0);
  tree->addMesh(mesh);

  center = 0.8125000000000000, 0.8125000000000000, 0.9375000000000000;
  halflengths = 0.0625000000000000, 0.0625000000000000, 0.0625000000000000;
  query::FindVoxelContent::BoundaryInclusion includeBounds =
      query::FindVoxelContent::INCLUDE_BOUNDARY;
  query::FindVoxelContent find(center, halflengths, includeBounds);

  //ExportSpacetree exporter("testSearchContentEdges-ReproduceBug");
  //precicePrint("-----------------------------------------------------------");
  //exporter.doExport(*tree);
  //precicePrint("+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++");
  //assertion(false);

  tree->searchContent(find);

  validateEquals(find.content().size(), 1);
}

void SpacetreeTestScenarios:: performTestSearchContentEdges
(
  int                     testDim,
  bool                    positive,
  const utils::DynVector& offset )
{
  preciceTrace3 ( "performTestSearchContentEdges()", testDim, positive, offset );
  int min = Spacetree::minElementsToRefineCell;
  Spacetree::minElementsToRefineCell = 1;
  int dim = offset.size();
  assertion ( not tarch::la::oneGreater(offset, utils::DynVector(dim,1.0)) );
  assertion ( tarch::la::allGreater(offset, utils::DynVector(dim,-1.0)) );
  bool flipNormals = false;
  mesh::PtrMesh mesh(new mesh::Mesh("TestMesh", dim, flipNormals));
  utils::DynVector coords0(offset);
  utils::DynVector coords1(offset);
  mesh::Vertex& v0 = mesh->createVertex(coords0);
  mesh::Vertex& v1 = mesh->createVertex(coords1);
  mesh::Vertex& v2 = mesh->createVertex(coords1);
  mesh::Vertex& v3 = mesh->createVertex(coords0);
  mesh->createEdge(v0, v1);
  mesh->createEdge(v2, v3);

  utils::DynVector center(dim, 0.0);
  utils::DynVector halflengths(dim, 1.0);
  query::FindVoxelContent::BoundaryInclusion includeBounds =
      query::FindVoxelContent::INCLUDE_BOUNDARY;
  query::FindVoxelContent::BoundaryInclusion excludeBounds =
      query::FindVoxelContent::EXCLUDE_BOUNDARY;
  query::FindVoxelContent findIncluded ( center, halflengths, includeBounds );
  query::FindVoxelContent findExcluded ( center, halflengths, excludeBounds );

  utils::DynVector treeoffset(dim, 0.0);
  utils::DynVector treeHalflengths(dim, 2.0);
  std::vector<double> refinementLimits;
  refinementLimits += 2.0, 1.0, 0.5, 0.25;
  std::vector<PtrSpacetree> treesInc;
  std::vector<PtrSpacetree> treesExc;
  for ( double limit : refinementLimits ) {
    treesInc += _factory.createSpacetree(treeoffset, treeHalflengths, limit);
    treesInc.back()->addMesh(mesh);
    treesExc += _factory.createSpacetree(treeoffset, treeHalflengths, limit);
    treesExc.back()->addMesh(mesh);
  }

  assertion(testDim >= 0);
  assertion(testDim < dim);

  double sign = positive ? 1.0 : -1.0;
  size_t size = 0;
  mesh::Merge merge;

  // Outside
  coords0[testDim] = sign * 1.5;
  coords1[testDim] = sign * 2.0;
  v0.setCoords(coords0);
  v1.setCoords(coords1);
  v2.setCoords(coords1);
  v3.setCoords(coords0);
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug("  outside i= " << i);
    preciceDebug("  including boundaries");
    treesInc[i]->searchContent(findIncluded);
    preciceDebug("  excluding boundaries");
    treesExc[i]->searchContent(findExcluded);
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent(findExcluded);
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().edges().size();
    validateEquals(size, 0);
    size = findExcluded.content().edges().size();
    validateEquals(size, 0);
  }

  // Outside eps
  coords0[testDim] = sign * (1.0 + 10.0 * tarch::la::NUMERICAL_ZERO_DIFFERENCE);
  coords1[testDim] = sign * 2.0;
  v0.setCoords ( coords0 );
  v1.setCoords ( coords1 );
  v2.setCoords(coords1);
  v3.setCoords(coords0);
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug ( "outside i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().edges().size();
    validateEquals ( size, 0 );
    size = findExcluded.content().edges().size();
    validateEquals ( size, 0 );
  }

  // Outside touching
  coords0[testDim] = sign * 1.0;
  coords1[testDim] = sign * 2.0;
  v0.setCoords ( coords0 );
  v1.setCoords ( coords1 );
  v2.setCoords(coords1);
  v3.setCoords(coords0);
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug ( "outside i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().edges().size();
    validateEquals ( size, 2 );
    findIncluded.clear ();
    size = findExcluded.content().edges().size();
    validateEquals ( size, 0 );
  }

  // Outside touching eps
  coords0[testDim] = sign * (1.0 + tarch::la::NUMERICAL_ZERO_DIFFERENCE);
  coords1[testDim] = sign * 2.0;
  v0.setCoords ( coords0 );
  v1.setCoords ( coords1 );
  v2.setCoords(coords1);
  v3.setCoords(coords0);
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug ( "outside i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().edges().size();
    validateEquals ( size, 2 );
    findIncluded.clear ();
    size = findExcluded.content().edges().size();
    validateEquals ( size, 0 );
  }

  // Intersecting
  coords0[testDim] = sign * 0.5;
  coords1[testDim] = sign * 1.5;
  v0.setCoords ( coords0 );
  v1.setCoords ( coords1 );
  v2.setCoords(coords1);
  v3.setCoords(coords0);
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug ( "outside i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().edges().size();
    validateEquals ( size, 2 );
    findIncluded.clear ();
    size = findExcluded.content().edges().size();
    validateEquals ( size, 2 );
    findExcluded.clear ();
  }

  // Intersecting eps
  coords0[testDim] = sign * (1.0 - 10.0 * tarch::la::NUMERICAL_ZERO_DIFFERENCE);
  coords1[testDim] = sign * 2.0;
  v0.setCoords ( coords0 );
  v1.setCoords ( coords1 );
  v2.setCoords(coords1);
  v3.setCoords(coords0);
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug ( "outside i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().edges().size();
    validateEquals ( size, 2 );
    findIncluded.clear ();
    size = findExcluded.content().edges().size();
    validateEquals ( size, 2 );
    findExcluded.clear ();
  }

  // Inside
  coords0[testDim] = sign * 0.3;
  coords1[testDim] = sign * 0.7;
  v0.setCoords ( coords0 );
  v1.setCoords ( coords1 );
  v2.setCoords(coords1);
  v3.setCoords(coords0);
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug ( "outside i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().edges().size();
    validateEquals ( size, 2 );
    findIncluded.clear ();
    size = findExcluded.content().edges().size();
    validateEquals ( size, 2 );
    findExcluded.clear ();
  }
  Spacetree::minElementsToRefineCell = min;
}

void SpacetreeTestScenarios:: testSearchContentTriangles()
{
  preciceTrace ( "testSearchContentTriangles()" );
  int dim = 3;
  for ( int testDim=0; testDim < dim; testDim++ ) {
    bool positiveDirection = true;
    bool negativeDirection = false;
    for ( int secondDim=0; secondDim < dim; secondDim++ ) {
      if ( secondDim == testDim ) {
        continue;
      }
      performTestSearchContentTriangles ( testDim, secondDim, positiveDirection );
      performTestSearchContentTriangles ( testDim, secondDim, negativeDirection );
    }
  }
}

void SpacetreeTestScenarios:: performTestSearchContentTriangles
(
  int  testDim,
  int  secondDimension,
  bool positive )
{
  preciceTrace2 ( "performTestSearchContentTriangles()", testDim, positive );
  int min = Spacetree::minElementsToRefineCell;
  Spacetree::minElementsToRefineCell = 1;
  int dim = 3;
  assertion ( testDim != secondDimension );
  bool flipNormals = false;
  mesh::PtrMesh mesh(new mesh::Mesh("TestMesh", dim, flipNormals));
  utils::Vector3D coords0(0.0);
  utils::Vector3D coords1(0.0);
  utils::Vector3D coords2(0.0);
  mesh::Vertex& v0 = mesh->createVertex(coords0);
  mesh::Vertex& v1 = mesh->createVertex(coords1);
  mesh::Vertex& v2 = mesh->createVertex(coords2);
  mesh::Edge& e0 = mesh->createEdge(v0, v1);
  mesh::Edge& e1 = mesh->createEdge(v1, v2);
  mesh::Edge& e2 = mesh->createEdge(v2, v0);
  mesh->createTriangle(e0, e1, e2);

  utils::DynVector center ( dim, 0.0 );
  utils::DynVector halflengths ( dim, 1.0 );
  query::FindVoxelContent::BoundaryInclusion includeBounds =
      query::FindVoxelContent::INCLUDE_BOUNDARY;
  query::FindVoxelContent::BoundaryInclusion excludeBounds =
      query::FindVoxelContent::EXCLUDE_BOUNDARY;
  query::FindVoxelContent findIncluded ( center, halflengths, includeBounds );
  query::FindVoxelContent findExcluded ( center, halflengths, excludeBounds );

  utils::DynVector treeoffset(dim, 0.0);
  utils::DynVector treeHalflengths(dim, 5.0);
  std::vector<double> refinementLimits;
  refinementLimits += 5.0, 2.5, 1.25, 0.625;
  std::vector<PtrSpacetree> treesInc;
  std::vector<PtrSpacetree> treesExc;
  for ( double limit : refinementLimits ) {
    treesInc += _factory.createSpacetree(treeoffset, treeHalflengths, limit);
    treesInc.back()->addMesh(mesh);
    treesExc += _factory.createSpacetree(treeoffset, treeHalflengths, limit);
    treesExc.back()->addMesh(mesh);
  }

  assertion ( testDim >= 0 );
  assertion ( testDim < dim );

  double sign = positive ? 1.0 : -1.0;
  size_t size = 0;
  mesh::Merge merge;

  int thirdDimension = -1;
  for ( int i=0; i < 3; i++ ) {
    if ( (i != testDim) && (i != secondDimension) ) {
      thirdDimension = i;
      break;
    }
  }

  // Outside
  coords0[testDim] = sign * 2.0;
  coords1[testDim] = sign * 3.0;
  coords2[testDim] = sign * 2.5;
  coords2[secondDimension] = sign * 0.5;
  v0.setCoords ( coords0 );
  v1.setCoords ( coords1 );
  v2.setCoords ( coords2 );
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug ( "i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().triangles().size();
    validateEquals ( size, 0 );
    size = findExcluded.content().triangles().size();
    validateEquals ( size, 0 );
  }

  // Outside eps vertex
  coords0[testDim] = sign * (1.0 + 10.0 * tarch::la::NUMERICAL_ZERO_DIFFERENCE);
  coords1[testDim] = sign * 2.0;
  coords2[testDim] = sign * 1.5;
  coords2[secondDimension] = sign * 0.5;
  v0.setCoords ( coords0 );
  v1.setCoords ( coords1 );
  v2.setCoords ( coords2 );
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug ( "i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().triangles().size();
    validateEquals ( size, 0 );
    size = findExcluded.content().triangles().size();
    validateEquals ( size, 0 );
  }

  // Outside eps edge
  coords0[testDim] = sign * (1.0 + 10.0 * tarch::la::NUMERICAL_ZERO_DIFFERENCE);
  coords1[testDim] = sign * 2.0;
  coords2[testDim] = sign * (1.0 + 10.0 * tarch::la::NUMERICAL_ZERO_DIFFERENCE);
  coords2[secondDimension] = sign * 0.5;
  v0.setCoords ( coords0 );
  v1.setCoords ( coords1 );
  v2.setCoords ( coords2 );
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug ( "i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().triangles().size();
    validateEquals ( size, 0 );
    size = findExcluded.content().triangles().size();
    validateEquals ( size, 0 );
  }

  // Touching eps vertex
  coords0[testDim] = sign * (1.0 + tarch::la::NUMERICAL_ZERO_DIFFERENCE);
  coords1[testDim] = sign * 2.0;
  coords2[testDim] = sign * 1.5;
  coords2[secondDimension] = sign * 0.5;
  v0.setCoords ( coords0 );
  v1.setCoords ( coords1 );
  v2.setCoords ( coords2 );
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug ( "i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().triangles().size();
    validateEquals ( size, 1 );
    findIncluded.clear ();
    size = findExcluded.content().triangles().size();
    validateEquals ( size, 0 );
  }

  // Touching eps edge
  coords0[testDim] = sign * (1.0 + tarch::la::NUMERICAL_ZERO_DIFFERENCE);
  coords1[testDim] = sign * 2.0;
  coords2[testDim] = sign * (1.0 + tarch::la::NUMERICAL_ZERO_DIFFERENCE);
  coords2[secondDimension] = sign * 0.5;
  v0.setCoords ( coords0 );
  v1.setCoords ( coords1 );
  v2.setCoords ( coords2 );
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug ( "i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().triangles().size();
    validateEquals ( size, 1 );
    findIncluded.clear ();
    size = findExcluded.content().triangles().size();
    validateEquals ( size, 0 );
  }

  // Touching vertex
  coords0[testDim] = sign * 1.0;
  coords1[testDim] = sign * 2.0;
  coords2[testDim] = sign * 1.5;
  coords2[secondDimension] = sign * 0.5;
  v0.setCoords ( coords0 );
  v1.setCoords ( coords1 );
  v2.setCoords ( coords2 );
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug ( "i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().triangles().size();
    validateEquals ( size, 1 );
    findIncluded.clear ();
    size = findExcluded.content().triangles().size();
    validateEquals ( size, 0 );
  }

  // Touching edge
  coords0[testDim] = sign * 1.0;
  coords1[testDim] = sign * 2.0;
  coords2[testDim] = sign * 1.0;
  coords2[secondDimension] = sign * 0.5;
  v0.setCoords ( coords0 );
  v1.setCoords ( coords1 );
  v2.setCoords ( coords2 );
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug ( "i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().triangles().size();
    validateEquals ( size, 1 );
    findIncluded.clear ();
    size = findExcluded.content().triangles().size();
    validateEquals ( size, 0 );
  }

  // Intersecting eps vertex
  coords0[testDim] = sign * (1.0 - 10.0 * tarch::la::NUMERICAL_ZERO_DIFFERENCE);
  coords1[testDim] = sign * 2.0;
  coords2[testDim] = sign * 1.5;
  coords2[secondDimension] = sign * 0.5;
  v0.setCoords ( coords0 );
  v1.setCoords ( coords1 );
  v2.setCoords ( coords2 );
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug ( "i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().triangles().size();
    validateEquals ( size, 1 );
    findIncluded.clear ();
    size = findExcluded.content().triangles().size();
    validateEquals ( size, 1 );
    findExcluded.clear ();
  }

  // Intersecting eps edge
  coords0[testDim] = sign * (1.0 - 10.0 * tarch::la::NUMERICAL_ZERO_DIFFERENCE);
  coords1[testDim] = sign * 2.0;
  coords2[testDim] = sign * (1.0 - 10.0 * tarch::la::NUMERICAL_ZERO_DIFFERENCE);
  coords2[secondDimension] = sign * 0.5;
  v0.setCoords ( coords0 );
  v1.setCoords ( coords1 );
  v2.setCoords ( coords2 );
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug ( "i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().triangles().size();
    validateEquals ( size, 1 );
    findIncluded.clear ();
    size = findExcluded.content().triangles().size();
    validateEquals ( size, 1 );
    findExcluded.clear ();
  }

  // Intersecting vertex
  coords0[testDim] = sign * 0.8;
  coords1[testDim] = sign * 2.0;
  coords2[testDim] = sign * 1.5;
  coords2[secondDimension] = sign * 0.5;
  v0.setCoords ( coords0 );
  v1.setCoords ( coords1 );
  v2.setCoords ( coords2 );
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug ( "i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().triangles().size();
    validateEquals ( size, 1 );
    findIncluded.clear ();
    size = findExcluded.content().triangles().size();
    validateEquals ( size, 1 );
    findExcluded.clear ();
  }

  // Intersecting edge
  coords0[testDim] = sign * 0.8;
  coords1[testDim] = sign * 2.0;
  coords2[testDim] = sign * 0.8;
  coords2[secondDimension] = sign * 0.5;
  v0.setCoords ( coords0 );
  v1.setCoords ( coords1 );
  v2.setCoords ( coords2 );
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug ( "i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().triangles().size();
    validateEquals ( size, 1 );
    findIncluded.clear ();
    size = findExcluded.content().triangles().size();
    validateEquals ( size, 1 );
    findExcluded.clear ();
  }

  // Contained
  coords0[testDim] = sign * 0.3;
  coords1[testDim] = sign * 0.8;
  coords2[testDim] = sign * 0.5;
  coords2[secondDimension] = sign * 0.5;
  v0.setCoords ( coords0 );
  v1.setCoords ( coords1 );
  v2.setCoords ( coords2 );
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug ( "i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().triangles().size();
    validateEquals ( size, 1 );
    findIncluded.clear ();
    size = findExcluded.content().triangles().size();
    validateEquals ( size, 1 );
    findExcluded.clear ();
  }

  // Contained filling
  coords0[testDim] = sign * -0.9;
  coords1[testDim] = sign * 0.9;
  coords2[testDim] = 0.0;
  coords2[secondDimension] = sign * 0.9;
  v0.setCoords ( coords0 );
  v1.setCoords ( coords1 );
  v2.setCoords ( coords2 );
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug ( "i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().triangles().size();
    validateEquals ( size, 1 );
    findIncluded.clear ();
    size = findExcluded.content().triangles().size();
    validateEquals ( size, 1 );
    findExcluded.clear ();
  }

  // Contained cutting
  coords0[testDim] = sign * -1.5;
  coords1[testDim] = sign * 1.5;
  coords2[testDim] = 0.0;
  coords2[secondDimension] = sign * 1.5;
  v0.setCoords ( coords0 );
  v1.setCoords ( coords1 );
  v2.setCoords ( coords2 );
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug ( "i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().triangles().size();
    validateEquals ( size, 1 );
    findIncluded.clear ();
    size = findExcluded.content().triangles().size();
    validateEquals ( size, 1 );
    findExcluded.clear ();
  }

  // Contained cutting wide1
  coords0[testDim] = sign * -5.0;
  coords1[testDim] = sign * 5.0;
  coords2[testDim] = 0.0;
  coords2[secondDimension] = sign * 5.0;
  v0.setCoords ( coords0 );
  v1.setCoords ( coords1 );
  v2.setCoords ( coords2 );
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug ( "i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().triangles().size();
    validateEquals ( size, 1 );
    findIncluded.clear ();
    size = findExcluded.content().triangles().size();
    validateEquals ( size, 1 );
    findExcluded.clear ();
  }

  // Touching contained fully
  coords0[testDim] = sign * 0.3;
  coords0[secondDimension] = sign * 1.0;
  coords1[testDim] = sign * 0.8;
  coords1[secondDimension] = sign * 1.0;
  coords2[testDim] = 0.5;
  coords2[secondDimension] = sign * 1.0;
  coords2[thirdDimension] = sign * 0.2;
  v0.setCoords ( coords0 );
  v1.setCoords ( coords1 );
  v2.setCoords ( coords2 );
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug ( "i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().triangles().size();
    validateEquals ( size, 1 );
    findIncluded.clear ();
    size = findExcluded.content().triangles().size();
    validateEquals ( size, 0 );
    findExcluded.clear ();
  }

  // Touching fully
  coords0[testDim] = sign * -1.5;
  coords0[secondDimension] = sign * 1.0;
  coords1[testDim] = sign * 1.5;
  coords1[secondDimension] = sign * 1.0;
  coords2[testDim] = 0.0;
  coords2[secondDimension] = sign * 1.0;
  coords2[thirdDimension] = sign * 1.5;
  v0.setCoords ( coords0 );
  v1.setCoords ( coords1 );
  v2.setCoords ( coords2 );
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug ( "i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().triangles().size();
    validateEquals ( size, 1 );
    findIncluded.clear ();
    size = findExcluded.content().triangles().size();
    validateEquals ( size, 0 );
    findExcluded.clear ();
  }

  // Touching fully wide
  coords0[testDim] = sign * -5.0;
  coords0[secondDimension] = sign * 1.0;
  coords1[testDim] = sign * 5.0;
  coords1[secondDimension] = sign * 1.0;
  coords2[testDim] = 0.0;
  coords2[secondDimension] = sign * 1.0;
  coords2[thirdDimension] = sign * 5.0;
  v0.setCoords ( coords0 );
  v1.setCoords ( coords1 );
  v2.setCoords ( coords2 );
  mesh->computeState();
  mesh->notifyListeners();
  for ( size_t i=0; i < treesInc.size(); i++ ) {
    preciceDebug ( "i= " << i );
    treesInc[i]->searchContent ( findIncluded );
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findIncluded.content() = merge(findIncluded.content());
    treesExc[i]->searchContent ( findExcluded );
    merge.content().clear();
    findExcluded.content() = merge(findExcluded.content());
    size = findIncluded.content().triangles().size();
    validateEquals ( size, 1 );
    findIncluded.clear ();
    size = findExcluded.content().triangles().size();
    validateEquals ( size, 0 );
    findExcluded.clear ();
  }
  Spacetree::minElementsToRefineCell = min;
}

void SpacetreeTestScenarios:: testVoxelPosition()
{
  preciceTrace ( "testVoxelPosition()" );
  int dim = 2;
  using utils::Vector2D;
  // geometry is a cuboid with offset (-3, -4) and sidelength 8
  utils::DynVector offset ( Vector2D(-3.0, -4.0) );
  utils::DynVector cuboidSidelength ( dim, 8.0 );
  bool flipNormals = true;
  mesh::PtrMesh mesh(new mesh::Mesh("test-cuboid", dim, flipNormals));
  geometry::Cuboid(offset, 0.01, cuboidSidelength).create(*mesh);

  io::ExportVTK exportVTK(true);
  exportVTK.doExport ( _testName + "-testVoxelPosition-cuboid", *mesh );

  //RegularSpacetree spacetree ( Vector(0.0, 0.0), 6.0, 0.05 );
  utils::DynVector treeOffset(dim, 0.0);
  utils::DynVector treeH(dim, 6.0);
  PtrSpacetree spacetree = _factory.createSpacetree ( treeOffset, treeH, 2.0 );
  spacetree->addMesh(mesh);
  spacetree->initialize();

  // voxels
  Vector2D voxelCenter1 (3.5, 1.5);
  Vector2D voxelHalflengths1 (1.5, 1.5);
  Vector2D voxelCenter2 (-3.0, 0.0);
  Vector2D voxelHalflengths2 (1.0, 1.0);
  Vector2D voxelCenter3 (-4.5, -4.0);
  Vector2D voxelHalflengths3 (0.5, 1.0);
  Vector2D voxelCenter4 (5.5, -4.5);
  Vector2D voxelHalflengths4 (0.5, 0.5);
  Vector2D voxelCenter5 (2.0, -2.0);
  Vector2D voxelHalflengths5 (1.0, 1.0);
  Vector2D voxelCenter6 (1.0, 0.0);
  Vector2D voxelHalflengths6 (1.0, 1.0);
  Vector2D voxelCenter7 (-2.6, 2.5);
  Vector2D voxelHalflengths7 (0.5, 0.5);
  Vector2D voxelCenter8 (-3.4, -2.5);
  Vector2D voxelHalflengths8 (0.5, 0.5);
  Vector2D voxelCenter9 (-4.5, 1.5);
  Vector2D voxelHalflengths9 (1.5, 1.5);

  query::FindVoxelContent voxel1 (
    voxelCenter1, voxelHalflengths1, query::FindVoxelContent::INCLUDE_BOUNDARY );
  int result1 = spacetree->searchContent ( voxel1 );

  query::FindVoxelContent voxel2 (
    voxelCenter2, voxelHalflengths2, query::FindVoxelContent::INCLUDE_BOUNDARY );
  int result2 = spacetree->searchContent ( voxel2 );

  query::FindVoxelContent voxel3 (
    voxelCenter3, voxelHalflengths3, query::FindVoxelContent::INCLUDE_BOUNDARY );
  int result3 = spacetree->searchContent ( voxel3 );

  query::FindVoxelContent voxel4 (
    voxelCenter4, voxelHalflengths4, query::FindVoxelContent::INCLUDE_BOUNDARY );
  int result4 = spacetree->searchContent ( voxel4 );

  query::FindVoxelContent voxel5 (
    voxelCenter5, voxelHalflengths5, query::FindVoxelContent::INCLUDE_BOUNDARY );
  int result5 = spacetree->searchContent ( voxel5 );

  query::FindVoxelContent voxel6 (
    voxelCenter6, voxelHalflengths6, query::FindVoxelContent::INCLUDE_BOUNDARY );
  int result6 = spacetree->searchContent ( voxel6 );

  query::FindVoxelContent voxel7 (
    voxelCenter7, voxelHalflengths7, query::FindVoxelContent::INCLUDE_BOUNDARY );
  int result7 = spacetree->searchContent ( voxel7 );

  query::FindVoxelContent voxel8 (
    voxelCenter8, voxelHalflengths8, query::FindVoxelContent::INCLUDE_BOUNDARY );
  int result8 = spacetree->searchContent ( voxel8 );

  query::FindVoxelContent voxel9 (
    voxelCenter9, voxelHalflengths9, query::FindVoxelContent::INCLUDE_BOUNDARY );
  int result9 = spacetree->searchContent ( voxel9 );

  ExportSpacetree export1 ( _testName + "testVoxelPosition-spacetree" );
  export1.doExport ( *spacetree );

//  query::ExportVTKVoxelQueries exportVoxels;
//  exportVoxels.addQuery ( voxelCenter1, voxelHalflengths1, voxel1.content().size() );
//  exportVoxels.addQuery ( voxelCenter2, voxelHalflengths2, voxel2.content().size() );
//  exportVoxels.addQuery ( voxelCenter3, voxelHalflengths3, voxel3.content().size() );
//  exportVoxels.addQuery ( voxelCenter4, voxelHalflengths4, voxel4.content().size() );
//  exportVoxels.addQuery ( voxelCenter5, voxelHalflengths5, voxel5.content().size() );
//  exportVoxels.addQuery ( voxelCenter6, voxelHalflengths6, voxel6.content().size() );
//  exportVoxels.addQuery ( voxelCenter7, voxelHalflengths7, voxel7.content().size() );
//  exportVoxels.addQuery ( voxelCenter8, voxelHalflengths8, voxel8.content().size() );
//  exportVoxels.addQuery ( voxelCenter9, voxelHalflengths9, voxel9.content().size() );
//
//  exportVoxels.exportQueries ("SpacetreeTest-testVoxelPosition-voxelqueries");

  validateEquals ( result1, Spacetree::positionOnGeometry() );
  validateEquals ( result2, Spacetree::positionOnGeometry() );
  validateEquals ( result3, Spacetree::positionInsideOfGeometry() );
  validateEquals ( result4, Spacetree::positionOnGeometry() );
  validateEquals ( result5, Spacetree::positionOutsideOfGeometry() );
  validateEquals ( result6, Spacetree::positionOutsideOfGeometry() );
  validateEquals ( result7, Spacetree::positionOnGeometry() );
  validateEquals ( result8, Spacetree::positionOnGeometry() );
  validateEquals ( result9, Spacetree::positionOnGeometry() );

//  // closest element from center of voxels to geometry
//  query::ExportVTKNeighbors exportNeighborsSpacetree;
//  Vector2D point1 (3.5, 1.5);
//  Vector2D point2 (-3.0, 0.0);
//  Vector2D point3 (-4.5, -4.0);
//  Vector2D point4 (5.5, -4.5);
//  Vector2D point5 (2.0, -2.0);
//  Vector2D point6 (1.0, 0.0);
//  Vector2D point7 (-2.6, 2.5);
//  Vector2D point8 (-3.4, -2.5);
//  Vector2D point9 (-4.5, 1.5);
//
//  query::FindClosest find1 ( point1 );
//  query::FindClosest find2 ( point2 );
//  query::FindClosest find3 ( point3 );
//  query::FindClosest find4 ( point4 );
//  query::FindClosest find5 ( point5 );
//  query::FindClosest find6 ( point6 );
//  query::FindClosest find7 ( point7 );
//  query::FindClosest find8 ( point8 );
//  query::FindClosest find9 ( point9 );
//
//  spacetree->searchDistance ( find1 );
//  spacetree->searchDistance ( find2 );
//  spacetree->searchDistance ( find3 );
//  spacetree->searchDistance ( find4 );
//  spacetree->searchDistance ( find5 );
//  spacetree->searchDistance ( find6 );
//  spacetree->searchDistance ( find7 );
//  spacetree->searchDistance ( find8 );
//  spacetree->searchDistance ( find9 );

//  exportNeighborsspacetree->addNeighbors ( point1, find1.getClosest() );
//  exportNeighborsspacetree->addNeighbors ( point2, find2.getClosest() );
//  exportNeighborsspacetree->addNeighbors ( point3, find3.getClosest() );
//  exportNeighborsspacetree->addNeighbors ( point4, find4.getClosest() );
//  exportNeighborsspacetree->addNeighbors ( point5, find5.getClosest() );
//  exportNeighborsspacetree->addNeighbors ( point6, find6.getClosest() );
//  exportNeighborsspacetree->addNeighbors ( point7, find7.getClosest() );
//  exportNeighborsspacetree->addNeighbors ( point8, find8.getClosest() );
//  exportNeighborsspacetree->addNeighbors ( point9, find9.getClosest() );

  ExportSpacetree export2 (_testName + "-testVoxelPosition-spacetree2");
  export2.doExport ( *spacetree );
  //exportNeighborsspacetree->exportNeighbors (
  //    "SpacetreeTest_testVoxelPosition-testVoxelPosition-withspacetree" );
}

void SpacetreeTestScenarios:: testSplittingVoxels()
{
  preciceTrace ( "testSplittingVoxels()" );
  int dim = 2;
  using utils::Vector2D;
  utils::DynVector cuboidOffset ( Vector2D(-5.0, -4.5) );
  utils::DynVector cuboidSidelength ( dim, 8.0 );
  bool flipNormals = true;
  mesh::PtrMesh mesh(new mesh::Mesh("test-cuboid", dim, flipNormals));
  geometry::Cuboid(cuboidOffset, 0.01, cuboidSidelength).create(*mesh);

  io::ExportVTK exportVTK(true);
  exportVTK.doExport ( _testName + "testSplittingVoxels-cuboid.vtk", *mesh );

  utils::DynVector treeOffset(dim, 0.0);
  utils::DynVector treeH(dim, 6.0);
  PtrSpacetree spacetree = _factory.createSpacetree(treeOffset, treeH, 0.05);
  spacetree->addMesh(mesh);
  spacetree->initialize();

  // voxels
  Vector2D voxelCenter1 (-3, 3.5);
  Vector2D voxelHalflengths1 (1, 1.5);
  Vector2D voxelCenter2 (2.5, 4.5);
  Vector2D voxelHalflengths2 (1.5, 0.5);
  Vector2D voxelCenter3 (-4.5, -0.5);
  Vector2D voxelHalflengths3 (0.5, 1.5);
  Vector2D voxelCenter4 (0.0, 0.5);
  Vector2D voxelHalflengths4 (2.0, 1.5);
  Vector2D voxelCenter5 (4.5, 0.0);
  Vector2D voxelHalflengths5 (1.5, 3.0);
  Vector2D voxelCenter6 (-3.75, -4.5);
  Vector2D voxelHalflengths6 (0.75, 1.5);
  Vector2D voxelCenter7 (-3, -4);
  Vector2D voxelHalflengths7 (1.0, 2.0);
  Vector2D voxelCenter8 (3.0, -4.5);
  Vector2D voxelHalflengths8 (2.0, 0.5);
  Vector2D voxelCenter9 (5.5, 5.0);
  Vector2D voxelHalflengths9 (0.5, 1.0);

  typedef Spacetree Spacetree;
  query::FindVoxelContent voxel1 (
      voxelCenter1, voxelHalflengths1, query::FindVoxelContent::INCLUDE_BOUNDARY );
  validateEquals ( spacetree->searchContent(voxel1), Spacetree::positionOnGeometry() );

  query::FindVoxelContent voxel2 (
      voxelCenter2, voxelHalflengths2, query::FindVoxelContent::INCLUDE_BOUNDARY );
  validateEquals ( spacetree->searchContent(voxel2), Spacetree::positionInsideOfGeometry() );

  query::FindVoxelContent voxel3 (
      voxelCenter3, voxelHalflengths3, query::FindVoxelContent::INCLUDE_BOUNDARY );
  validateEquals ( spacetree->searchContent(voxel3), Spacetree::positionOnGeometry() );

  query::FindVoxelContent voxel4 (
      voxelCenter4, voxelHalflengths4, query::FindVoxelContent::INCLUDE_BOUNDARY );
  validateEquals ( spacetree->searchContent(voxel4), Spacetree::positionOutsideOfGeometry() );

  query::FindVoxelContent voxel5 (
      voxelCenter5, voxelHalflengths5, query::FindVoxelContent::INCLUDE_BOUNDARY );
  validateEquals ( spacetree->searchContent(voxel5), Spacetree::positionOnGeometry() );

  query::FindVoxelContent voxel6 (
      voxelCenter6, voxelHalflengths6, query::FindVoxelContent::INCLUDE_BOUNDARY );
  validateEquals (
      spacetree->searchContent(voxel6), Spacetree::positionOnGeometry() );

  query::FindVoxelContent voxel7 (
      voxelCenter7, voxelHalflengths7, query::FindVoxelContent::INCLUDE_BOUNDARY );
  validateEquals ( spacetree->searchContent(voxel7), Spacetree::positionOnGeometry() );

  query::FindVoxelContent voxel8 (
      voxelCenter8, voxelHalflengths8, query::FindVoxelContent::INCLUDE_BOUNDARY );
  validateEquals ( spacetree->searchContent(voxel8), Spacetree::positionOnGeometry() );

  query::FindVoxelContent voxel9 (
      voxelCenter9, voxelHalflengths9, query::FindVoxelContent::INCLUDE_BOUNDARY );
  validateEquals ( spacetree->searchContent(voxel9), Spacetree::positionInsideOfGeometry() );
}

}}} // namespace precice, spacetree, tests

