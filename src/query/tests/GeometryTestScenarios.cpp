#include "GeometryTestScenarios.hpp"
#include "precice/Constants.hpp"

namespace precice {
namespace query {
namespace tests {

GeometryTestScenarios:: GeometryTestScenarios()
:
  _pointQueryScenario ( nullptr ),
  _positionQueryScenario ( nullptr ),
  _voxelQueryScenario ( nullptr )
{}

GeometryTestScenarios:: ~GeometryTestScenarios()
{
  if ( _pointQueryScenario != nullptr ){
    delete _pointQueryScenario;
    _pointQueryScenario = nullptr;
  }
  if ( _positionQueryScenario != nullptr ){
    delete _positionQueryScenario;
    _positionQueryScenario = nullptr;
  }
  if ( _voxelQueryScenario != nullptr ){
    delete _voxelQueryScenario;
    _voxelQueryScenario = nullptr;
  }
}

GeometryTestScenarios::PointQueryScenario:: PointQueryScenario
(
  int dim )
:
  mesh("PointQueryScenarioMesh", dim, false),
  queryCoords(),
  validDistances(),
  validDistanceVectors()
{}

GeometryTestScenarios::PositionQueryScenario:: PositionQueryScenario
(
  int dim )
:
  mesh("PositionQueryScenarioMesh", dim, false),
  queryCoords(),
  validPositions()
{}

GeometryTestScenarios::VoxelQueryScenario:: VoxelQueryScenario
(
  int dim )
:
  mesh("VoxelQueryScenarioMesh", dim, false),
  queryCenters(),
  queryHalflengths(),
  includeBoundaries(),
  validPositions()
{}

const GeometryTestScenarios::PointQueryScenario&
GeometryTestScenarios:: pointQueryScenario
(
  int dim )
{
  using utils::DynVector;
  if ( _pointQueryScenario == nullptr ){
    _pointQueryScenario = new PointQueryScenario(dim);

    // Create query mesh (square in 2D, cube in 3D)
    mesh::Mesh& mesh = _pointQueryScenario->mesh;
    createMesh(mesh);

    // Create query positions and valid results
    for ( int testDim=0; testDim < dim; testDim++ ){
      DynVector coord ( dim, 0.0 );
      double distance = 0.0;
      DynVector vector ( dim, 0.0 );

      coord[testDim] = 0.1;
      distance = 0.9;
      vector[testDim] = 0.9;
      addToPointQueryScenario ( coord, distance, vector);

      coord[testDim] = 0.5;
      distance = 0.5;
      vector[testDim] = 0.5;
      addToPointQueryScenario ( coord, distance, vector );

      coord[testDim] = 0.9;
      distance = 0.1;
      vector[testDim] = 0.1;
      addToPointQueryScenario ( coord, distance, vector );
    }
  }
  return *_pointQueryScenario;
}

const GeometryTestScenarios::PositionQueryScenario &
GeometryTestScenarios:: positionQueryScenario
(
  int dim )
{
  using utils::DynVector;
  if ( _positionQueryScenario == nullptr ){
    _positionQueryScenario = new PositionQueryScenario(dim);

    // Create query mesh (square in 2D, cube in 3D)
    mesh::Mesh& mesh = _positionQueryScenario->mesh;
    createMesh ( mesh );

    // Create query positions and valid results
    std::list<DynVector>& coords = _positionQueryScenario->queryCoords;
    std::list<int>& positions = _positionQueryScenario->validPositions;

    for ( int testDim=0; testDim < dim; testDim++ ){
      DynVector coord ( dim, 0.0 );

      coord[testDim] = 0.1;
      coords.push_back ( coord );
      positions.push_back ( constants::positionInsideOfGeometry() );

      coord[testDim] = 0.5;
      coords.push_back ( coord );
      positions.push_back ( constants::positionInsideOfGeometry() );

      coord[testDim] = 0.9;
      coords.push_back ( coord );
      positions.push_back ( constants::positionInsideOfGeometry() );
    }
  }
  return *_positionQueryScenario;
}

const GeometryTestScenarios::VoxelQueryScenario &
GeometryTestScenarios:: voxelQueryScenario
(
  int dim )
{
  using utils::DynVector;
  if ( _voxelQueryScenario == nullptr ){
    _voxelQueryScenario = new VoxelQueryScenario(dim);

    // Create query mesh (square in 2D, cube in 3D)
    mesh::Mesh& mesh = _voxelQueryScenario->mesh;
    createMesh ( mesh );

    // Create query positions and valid results
    std::list<DynVector>& queryCenters = _voxelQueryScenario->queryCenters;
    std::list<DynVector>& queryHalflengths = _voxelQueryScenario->queryHalflengths;
    std::list<bool>& includeBoundaries = _voxelQueryScenario->includeBoundaries;
    std::list<int>& positions = _voxelQueryScenario->validPositions;
    for ( int testDim=0; testDim < dim; testDim++ ){
      DynVector center ( dim, 0.0 );
      DynVector halflengths ( dim, 0.1 );

      center[testDim] = 0.1;
      queryCenters.push_back ( center );
      queryHalflengths.push_back ( halflengths );
      includeBoundaries.push_back ( false );
      positions.push_back ( constants::positionInsideOfGeometry() );

      center[testDim] = 0.5;
      queryCenters.push_back ( center );
      queryHalflengths.push_back ( halflengths );
      includeBoundaries.push_back ( false );
      positions.push_back ( constants::positionInsideOfGeometry() );

      center[testDim] = 0.9;
      queryCenters.push_back ( center );
      queryHalflengths.push_back ( halflengths );
      includeBoundaries.push_back ( false );
      positions.push_back ( constants::positionInsideOfGeometry() );
    }
  }
  return *_voxelQueryScenario;
}

void GeometryTestScenarios:: createMesh
(
  mesh::Mesh& mesh )
{
  using namespace mesh;
  if (mesh.getDimensions() == 2){
    Vertex& v00 = mesh.createVertex(Eigen::Vector2d(-1.0, -1.0));
    Vertex& v10 = mesh.createVertex(Eigen::Vector2d( 1.0, -1.0));
    Vertex& v01 = mesh.createVertex(Eigen::Vector2d(-1.0,  1.0));
    Vertex& v11 = mesh.createVertex(Eigen::Vector2d( 1.0,  1.0));
    mesh.createEdge(v00, v10);
    mesh.createEdge(v10, v11);
    mesh.createEdge(v11, v01);
    mesh.createEdge(v01, v00);
  }
  else {
    assertion(mesh.getDimensions() == 3, mesh.getDimensions());
    Vertex& v000 = mesh.createVertex(Eigen::Vector3d(-1.0, -1.0, -1.0));
    Vertex& v001 = mesh.createVertex(Eigen::Vector3d(-1.0, -1.0,  1.0));
    Vertex& v010 = mesh.createVertex(Eigen::Vector3d(-1.0,  1.0, -1.0)); //
    Vertex& v011 = mesh.createVertex(Eigen::Vector3d(-1.0,  1.0,  1.0));
    Vertex& v100 = mesh.createVertex(Eigen::Vector3d( 1.0, -1.0, -1.0)); //
    Vertex& v101 = mesh.createVertex(Eigen::Vector3d( 1.0, -1.0,  1.0));
    Vertex& v110 = mesh.createVertex(Eigen::Vector3d( 1.0,  1.0, -1.0));
    Vertex& v111 = mesh.createVertex(Eigen::Vector3d( 1.0,  1.0,  1.0));

    Edge& e000to100 = mesh.createEdge(v000, v100);
    Edge& e010to110 = mesh.createEdge(v010, v110);
    Edge& e001to101 = mesh.createEdge(v001, v101);
    Edge& e011to111 = mesh.createEdge(v011, v111);

    Edge& e000to010 = mesh.createEdge(v000, v010);
    Edge& e100to110 = mesh.createEdge(v100, v110);
    Edge& e001to011 = mesh.createEdge(v001, v011);
    Edge& e101to111 = mesh.createEdge(v101, v111);

    Edge& e000to001 = mesh.createEdge(v000, v001);
    Edge& e100to101 = mesh.createEdge(v100, v101);
    Edge& e010to011 = mesh.createEdge(v010, v011);
    Edge& e110to111 = mesh.createEdge(v110, v111);
    Edge& e010to001 = mesh.createEdge(v010, v001);

    Edge& e100to111 = mesh.createEdge(v100, v111);
    Edge& e100to001 = mesh.createEdge(v100, v001);
    Edge& e110to011 = mesh.createEdge(v110, v011);
    Edge& e010to100 = mesh.createEdge(v010, v100); //
    Edge& e001to111 = mesh.createEdge(v001, v111);

    mesh.createTriangle(e000to001, e010to001, e000to010); // x = 0
    mesh.createTriangle(e010to011, e010to001, e001to011);
    mesh.createTriangle(e100to101, e100to111, e101to111); // x = 1
    mesh.createTriangle(e100to110, e110to111, e100to111);

    mesh.createTriangle(e000to100, e100to001, e000to001); // y = 0
    mesh.createTriangle(e100to101, e001to101, e100to001);
    mesh.createTriangle(e010to110, e010to011, e110to011); // y = 1
    mesh.createTriangle(e110to111, e110to011, e011to111);

    mesh.createTriangle(e000to100, e000to010, e010to100); // z = 0
    mesh.createTriangle(e010to100, e010to110, e100to110);
    mesh.createTriangle(e001to101, e101to111, e001to111); // z = 1
    mesh.createTriangle(e001to011, e001to111, e011to111);
  }
  mesh.computeState();
}

void GeometryTestScenarios:: addToPointQueryScenario
(
  const utils::DynVector& queryCoord,
  double                  validDistance,
  const utils::DynVector& validDistanceVector )
{
  assertion ( _pointQueryScenario != nullptr );
  _pointQueryScenario->queryCoords.push_back ( queryCoord );
  _pointQueryScenario->validDistances.push_back ( validDistance );
  _pointQueryScenario->validDistanceVectors.push_back ( validDistanceVector );
}

}}}// namespace precice, query, tests
