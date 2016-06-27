#include "FindVoxelContent.hpp"
#include "utils/GeometryComputations.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Mesh.hpp"
#include "tarch/la/VectorOperations.h"
#include "tarch/la/VectorVectorOperations.h"

namespace precice {
namespace query {

//#define PRECICE_OLD_QUERY

logging::Logger FindVoxelContent::_log("precice::query::FindVoxelContent");

//FindVoxelContent:: FindVoxelContent
//(
//  const utils::DynVector& voxelCenter,
//  const utils::DynVector& halflengths,
//  BoundaryInclusion       boundaryInclusion )
//:
//  _voxelCenter ( voxelCenter ),
//  _voxelHalflengths ( halflengths ),
//  _boundaryInclusion ( boundaryInclusion ),
//  _dimensions ( voxelCenter.size() ),
//  _content ()
//{
//  assertion ( voxelCenter.size() == halflengths.size(),
//               voxelCenter.size(), halflengths.size() );
//  assertion ( (_dimensions == 2) || (_dimensions == 3), _dimensions );
//}

//void FindVoxelContent:: setVoxelCenter
//(
//  const utils::DynVector& center )
//{
//  _voxelCenter = center;
//}

//void FindVoxelContent:: setVoxelHalflengths
//(
//  const utils::DynVector& halflengths )
//{
//  _voxelHalflengths = halflengths;
//}

const utils::DynVector& FindVoxelContent:: getVoxelCenter() const
{
  return _voxelCenter;
}

const utils::DynVector& FindVoxelContent:: getVoxelHalflengths() const
{
  return _voxelHalflengths;
}

FindVoxelContent::BoundaryInclusion FindVoxelContent:: getBoundaryInclusion() const
{
  return _boundaryInclusion;
}

mesh::Group& FindVoxelContent:: content()
{
  return _content;
}

void FindVoxelContent:: clear()
{
  _content.clear();
}

void FindVoxelContent:: checkVertex
(
  mesh::Vertex& vertex )
{
  using namespace tarch::la;
  assertion(vertex.getDimensions() == _dimensions, vertex.getDimensions(),
             _dimensions);
  utils::DynVector toVertex = vertex.getCoords();
  toVertex -= _voxelCenter;
  abs(toVertex, toVertex);
  utils::DynVector absHalflengths(_voxelHalflengths.size(), 0.0);
  abs(_voxelHalflengths, absHalflengths);
  if (_boundaryInclusion == INCLUDE_BOUNDARY){
    if (oneGreater(toVertex, absHalflengths)){
      return;
    }
  }
  else if (oneGreaterEquals(toVertex, absHalflengths)){
    return;
  }
//#   ifdef Debug
//    if (_voxelCenter.size() == 3)
//    if (norm2(_voxelCenter - utils::Vector3D(4.343209876543000, 4.0666666666666666, 4.106172839509999)) < 1e-10){
//      precicePrint("Add vertex " << vertex.getCoords());
//    }
//#   endif
  _content.add(vertex);
}

void FindVoxelContent:: checkEdge
(
  mesh::Edge& edge )
{
  preciceTrace("checkEdge()", edge.vertex(0).getCoords(), edge.vertex(1).getCoords());
  using namespace tarch::la;
  using utils::Vector3D;
  using utils::Vector2D;
  using utils::DynVector;
  assertion(edge.getDimensions() == _dimensions, edge.getDimensions(),
             _dimensions);
# ifdef PRECICE_OLD_QUERY
  DynVector toEdgeCenter = edge.getCenter();
  toEdgeCenter -= _voxelCenter;
  abs(toEdgeCenter, toEdgeCenter);
  // Test, if circumcircle (circumsphere) of edge lies completely outside
  if (oneGreater(toEdgeCenter - edge.getEnclosingRadius(), _voxelHalflengths)){
    preciceDebug("  edge circumcircle lies outside");
    return;
  }
  // Test, if edge center lies completely inside
  DynVector absVoxelHalflengths(_dimensions);
  if (allGreater(tarch::la::abs(_voxelHalflengths, absVoxelHalflengths), toEdgeCenter)){
    _content.add(edge);
    preciceDebug("  edge center lies inside");
    return;
  }

# else // PRECICE_OLD_QUERY
  // For excluding boundary case, add eps to achieve greater/smaller-equals
  double eps = 0.0;
  if ( _boundaryInclusion == EXCLUDE_BOUNDARY ){
    eps = 2.0 * NUMERICAL_ZERO_DIFFERENCE;
  }
# endif // not PRECICE_OLD_QUERY

  if ( _dimensions == 2 ){
#   ifndef PRECICE_OLD_QUERY
    // Test if edge intersects with rectangle, using seperating axis theorem
    Vector2D a = edge.vertex(0).getCoords();
    a -= _voxelCenter; // Move to origin
    Vector2D b = edge.vertex(1).getCoords();
    b -= _voxelCenter; // Move to origin
    for (int dim=0; dim < 2; dim++ ){
      // Test projection in "dim"
      double edgeMax = a[dim] > b[dim] ? a[dim] : b[dim];
      double edgeMin = a[dim] < b[dim] ? a[dim] : b[dim];
      double voxelMax = _voxelHalflengths[dim];
      double voxelMin = -1.0 * _voxelHalflengths[dim];
      if ( greater(edgeMin, voxelMax-eps) || smaller(edgeMax-eps, voxelMin) ){
        return;
      }
    }
    // Test projection in edge normal plane
    // Further remove suitable dimension to work in 1D
    int dim = indexMax(abs(edge.getNormal())); // 1D dimensions
    double projEdge = dot(edge.getNormal(), a) * edge.getNormal()[dim];
    tarch::la::Vector<4,double> projVoxel;
    Vector2D corner ( -_voxelHalflengths[0], -_voxelHalflengths[1] );
    projVoxel[0] = dot(edge.getNormal(), corner) * edge.getNormal()[dim];
    assignList(corner) = _voxelHalflengths[0], -_voxelHalflengths[1];
    projVoxel[1] = dot(edge.getNormal(), corner) * edge.getNormal()[dim];
    assignList(corner) = -_voxelHalflengths[0], _voxelHalflengths[1];
    projVoxel[2] = dot(edge.getNormal(), corner) * edge.getNormal()[dim];
    assignList(corner) = _voxelHalflengths[0], _voxelHalflengths[1];
    projVoxel[3] = dot(edge.getNormal(), corner) * edge.getNormal()[dim];
    double voxelMax = max(projVoxel);
    double voxelMin = min(projVoxel);
    if ( greater(projEdge, voxelMax-eps) || smaller(projEdge-eps, voxelMin) ){
      return;
    }
    // If all tests pass, the edge intersects the voxel
    _content.add (edge);

#   else // not PRECICE_OLD_QUERY
    // Test, if edge intersects with voxel
    Vector2D a = edge.vertex(0).getCoords();
    Vector2D b = edge.vertex(1).getCoords();
    // lower edge of voxel
    Vector2D voxelPointA, voxelPointB;
    voxelPointA = _voxelCenter;
    voxelPointA -= _voxelHalflengths;
    assignList(voxelPointB) = _voxelCenter[0] + _voxelHalflengths[0],
                              _voxelCenter[1] - _voxelHalflengths[1];
    if ( utils::GeometryComputations::segmentsIntersect(
         voxelPointA, voxelPointB, a, b, _boundaryInclusion == INCLUDE_BOUNDARY) )
    {
      _content.add (edge);
      preciceDebug ( "  lower edge of voxel intersects edge" );
      return;
    }
    // right edge
    voxelPointA = voxelPointB;
    voxelPointB = _voxelCenter;
    voxelPointB += _voxelHalflengths;
    if ( utils::GeometryComputations::segmentsIntersect(
         voxelPointA, voxelPointB, a, b, _boundaryInclusion == INCLUDE_BOUNDARY) )
    {
      _content.add (edge);
      preciceDebug ( "  right edge of voxel intersects edge" );
      return;
    }
    // top edge
    voxelPointA = voxelPointB;
    assignList(voxelPointB) = _voxelCenter[0] - _voxelHalflengths[0],
                              _voxelCenter[1] + _voxelHalflengths[1];
    if ( utils::GeometryComputations::segmentsIntersect(
         voxelPointA, voxelPointB, a, b, _boundaryInclusion == INCLUDE_BOUNDARY) )
    {
      _content.add (edge);
      preciceDebug ( "  top edge of voxel intersects edge" );
      return;
    }
    // left edge
    voxelPointA = voxelPointB;
    voxelPointB = _voxelCenter;
    voxelPointB -= _voxelHalflengths;
    if ( utils::GeometryComputations::segmentsIntersect(
         voxelPointA, voxelPointB, a, b, _boundaryInclusion == INCLUDE_BOUNDARY) )
    {
      _content.add (edge);
      preciceDebug ( "  left edge of voxel intersects edge" );
      return;
    }
#   endif // PRECIC_OLD_QUERY
  }
  else { // 3D
    assertion(_dimensions == 3, _dimensions);
#   ifndef PRECICE_OLD_QUERY
    // Test if edge intersects with rectangle, using seperating axis theorem
    Vector3D a = edge.vertex(0).getCoords();
    a -= _voxelCenter; // Move to origin
    Vector3D b = edge.vertex(1).getCoords();
    b -= _voxelCenter; // Move to origin
    for (int dim=0; dim < 3; dim++){
      // Test projection in "dim"
      double edgeMax = a[dim] > b[dim] ? a[dim] : b[dim];
      double edgeMin = a[dim] < b[dim] ? a[dim] : b[dim];
      double voxelMax = _voxelHalflengths[dim];
      double voxelMin = -1.0 * _voxelHalflengths[dim];
      if (greater(edgeMin, voxelMax-eps) || smaller(edgeMax-eps, voxelMin)){
        return;
      }
    }

    // Test projection in edge normal plane
    // Further remove suitable dimension to work in 1D
//    utils::Vector3D absNormal;
//    int dim = indexMax(abs(edge.getNormal(), absNormal)); // 1D dimensions
//    double projEdge = dot(edge.getNormal(), a) * edge.getNormal()[dim];
//    tarch::la::Vector<8,double> projVoxel;
//    utils::DynVector& h = _voxelHalflengths;
//    utils::Vector3D corner ( -h[0], -h[1], -h[2] );
//    projVoxel[0] = dot(edge.getNormal(), corner) * edge.getNormal()[dim];
//    assignList(corner) =  h[0], -h[1], -h[2];
//    projVoxel[1] = dot(edge.getNormal(), corner) * edge.getNormal()[dim];
//    assignList(corner) = -h[0],  h[1], -h[2];
//    projVoxel[2] = dot(edge.getNormal(), corner) * edge.getNormal()[dim];
//    assignList(corner) =  h[0],  h[1], -h[2];
//    projVoxel[3] = dot(edge.getNormal(), corner) * edge.getNormal()[dim];
//    assignList(corner) = -h[0], -h[1],  h[2];
//    projVoxel[4] = dot(edge.getNormal(), corner) * edge.getNormal()[dim];
//    assignList(corner) =  h[0], -h[1],  h[2];
//    projVoxel[5] = dot(edge.getNormal(), corner) * edge.getNormal()[dim];
//    assignList(corner) = -h[0],  h[1],  h[2];
//    projVoxel[6] = dot(edge.getNormal(), corner) * edge.getNormal()[dim];
//    assignList(corner) =  h[0],  h[1],  h[2];
//    projVoxel[7] = dot(edge.getNormal(), corner) * edge.getNormal()[dim];
//    double voxelMax = max(projVoxel);
//    double voxelMin = min(projVoxel);
//    if ( greater(projEdge, voxelMax-eps) || smaller(projEdge-eps, voxelMin) ){
//      return;
//    }
    Vector3D edgeVector = b;
    edgeVector -= a;
    // Separation axis is formed by cross product of unit vector in x-direction
    // and edge vector:
    Vector3D sepAxis(0.0, -edgeVector[2], edgeVector[1]); // Dim 0
    // Edge is projected onto separation axis by dot product:
    double proj = sepAxis[1]*a[1] + sepAxis[2]*a[2];
    double radius = _voxelHalflengths[1]*abs(sepAxis[1])
                    + _voxelHalflengths[2]*abs(sepAxis[2]);
    if (greater(proj, radius-eps) || smaller(proj-eps, -radius)){
      if (radius > 0.0) return;
    }
    assignList(sepAxis) = edgeVector[2], 0.0, -edgeVector[0]; // Dim 1
    proj = sepAxis[0]*a[0] + sepAxis[2]*a[2];
    radius = _voxelHalflengths[0]*abs(sepAxis[0])
             + _voxelHalflengths[2]*abs(sepAxis[2]);
    if (greater(proj, radius-eps) || smaller(proj-eps, -radius)){
      if (radius > 0.0) return;
    }
    assignList(sepAxis) = -edgeVector[1], edgeVector[0], 0.0; // Dim 2
    proj = sepAxis[0]*a[0] + sepAxis[1]*a[1];
    radius = _voxelHalflengths[0]*abs(sepAxis[0])
             + _voxelHalflengths[1]*abs(sepAxis[1]);
    if (greater(proj, radius-eps) || smaller(proj-eps, -radius)){
      if (radius > 0.0) return;
    }

//#   ifdef Debug
//    if (_voxelCenter.size() == 3)
//    if (norm2(_voxelCenter - Vector3D(4.343209876543000, 4.0666666666666666, 4.106172839509999)) < 1e-10){
//      precicePrint("Add edge from " << edge.vertex(0).getCoords() << " to " << edge.vertex(1).getCoords());
//    }
//#   endif

    // If all tests pass, the edge intersects the voxel
    _content.add(edge);
#   else // not PRECICE_OLD_QUERY
    typedef utils::GeometryComputations Geocomp;
    bool include = _boundaryInclusion == INCLUDE_BOUNDARY;

    Vector3D coords0;
    Vector3D coords1;
    coords0 = edge.vertex(0).getCoords();
    coords1 = edge.vertex(1).getCoords();
    Vector3D center;
    // Voxel face orthogonal to x (at x < center)
    assignList(center) = _voxelCenter[0] - _voxelHalflengths[0], _voxelCenter[1], _voxelCenter[2];
    if ( computeIntersection(center, _voxelHalflengths, 0, coords0, coords1, include) ) {
      _content.add ( edge );
      return;
    }
    // Voxel face orthogonal to x (at x > center)
    assignList(center) = _voxelCenter[0] + _voxelHalflengths[0], _voxelCenter[1], _voxelCenter[2];
    if ( computeIntersection(center, _voxelHalflengths, 0, coords0, coords1, include) ) {
      _content.add ( edge );
      return;
    }
    // Voxel face orthogonal to y (at y < center)
    assignList(center) = _voxelCenter[0], _voxelCenter[1] - _voxelHalflengths[1], _voxelCenter[2];
    if ( computeIntersection(center, _voxelHalflengths, 1, coords0, coords1, include) ) {
      _content.add ( edge );
      return;
    }
    // Voxel face orthogonal to y (at y > center)
    assignList(center) = _voxelCenter[0], _voxelCenter[1] + _voxelHalflengths[1], _voxelCenter[2];
    if ( computeIntersection(center, _voxelHalflengths, 1, coords0, coords1, include) ) {
      _content.add ( edge );
      return;
    }
    // Voxel face orthogonal to z (at z < center)
    assignList(center) = _voxelCenter[0], _voxelCenter[1], _voxelCenter[2] - _voxelHalflengths[2];
    if ( computeIntersection(center, _voxelHalflengths, 2, coords0, coords1, include) ) {
      _content.add ( edge );
      return;
    }
    // Voxel face orthogonal to z (at z > center)
    assignList(center) = _voxelCenter[0], _voxelCenter[1], _voxelCenter[2] + _voxelHalflengths[2];
    if ( computeIntersection(center, _voxelHalflengths, 2, coords0, coords1, include) ) {
      _content.add ( edge );
      return;
    }

    for ( int i=0; i < 2; i ++ ) {
      DynVector toEdge = _voxelCenter;
      toEdge -= edge.vertex(i).getCoords();
      abs ( toEdge, toEdge );
      bool touching = not oneGreater ( toEdge, _voxelHalflengths );
      touching &= not allGreater ( _voxelHalflengths, toEdge );
      DynVector edgeDirection = edge.vertex(1-i).getCoords();
      edgeDirection -= edge.vertex(i).getCoords();
      edgeDirection /= norm2 ( edgeDirection );
      edgeDirection *= 10.0 * NUMERICAL_ZERO_DIFFERENCE;
      toEdge = _voxelCenter;
      toEdge -= edge.vertex(i).getCoords();
      toEdge -= edgeDirection;
      abs ( toEdge, toEdge );
      if ( allGreater(_voxelHalflengths, toEdge) ) {
        _content.add ( edge );
        preciceDebug ( "  edge completely intersects voxel" );
        return;
      }
    }
#   endif // PRECICE_OLD_QUERY
  }
}

void FindVoxelContent:: checkTriangle
(
  mesh::Triangle& triangle )
{
  preciceTrace ( "checkTriangle()", triangle.getID(), triangle.getCenter() );
  assertion ( _dimensions == 3, _dimensions );
  using namespace tarch::la;
  using utils::Vector3D;
  using utils::Vector2D;

# ifndef PRECICE_OLD_QUERY
  // Use Separating axis theorem to find intersection

  // For excluding boundary case, add eps to achieve greater/smaller-equals
  double eps = 0.0;
  if ( _boundaryInclusion == EXCLUDE_BOUNDARY ){
    eps = 2.0 * NUMERICAL_ZERO_DIFFERENCE;
  }

  // Shift triangle vertex coordinates to have voxel at coordinate origin
  tarch::la::Vector<3,Vector3D> vertices;
  for (int i=0; i < 3; i++){
    vertices[i] = triangle.vertex(i).getCoords();
    vertices[i] -= _voxelCenter;
  }

  // 1. Tests:
  // Voxel sides against minimal triangle axis-aligned bounding box
  Vector3D coords; // Holds triangle coords of one dimension
  double triangleMin = 0.0;
  double triangleMax = 0.0;
  double voxelMax = 0.0;
  double voxelMin = 0.0;
  for (int dim=0; dim < 3; dim++ ){     // Test projection in "dim"
    for (int i=0; i < 3; i++){
      coords[i] = vertices[i][dim];
    }
    triangleMin = min(coords);
    voxelMax = _voxelHalflengths[dim];
    if ( greater(triangleMin, voxelMax-eps) ){
      preciceDebug ( "Found sa in test 1.1" );
      return;
    }
    triangleMax = max(coords);
    voxelMin = -1.0 * _voxelHalflengths[dim];
    if ( smaller(triangleMax-eps, voxelMin) ){
      preciceDebug ( "Found sa in test 1.2" );
      return;
    }
  }

  // 2. Tests:
  // Triangle plane intersecting with voxel. Projection of voxel vertices on
  // triangle normal.
  utils::Vector3D n = triangle.getNormal();
  int dim = indexMax(abs(n)); // 1D dimensions
  double projTri = dot(n, vertices[0]) * n[dim];
  tarch::la::Vector<8,double> projVoxel;
  utils::DynVector& h = _voxelHalflengths;
  utils::Vector3D corner ( -h[0], -h[1], -h[2] );
  projVoxel[0] = dot(n, corner) * n[dim];
  assignList(corner) =  h[0], -h[1], -h[2];
  projVoxel[1] = dot(n, corner) * n[dim];
  assignList(corner) = -h[0],  h[1], -h[2];
  projVoxel[2] = dot(n, corner) * n[dim];
  assignList(corner) =  h[0],  h[1], -h[2];
  projVoxel[3] = dot(n, corner) * n[dim];
  assignList(corner) = -h[0], -h[1],  h[2];
  projVoxel[4] = dot(n, corner) * n[dim];
  assignList(corner) =  h[0], -h[1],  h[2];
  projVoxel[5] = dot(n, corner) * n[dim];
  assignList(corner) = -h[0],  h[1],  h[2];
  projVoxel[6] = dot(n, corner) * n[dim];
  assignList(corner) =  h[0],  h[1],  h[2];
  projVoxel[7] = dot(n, corner) * n[dim];
  voxelMax = max(projVoxel);
  voxelMin = min(projVoxel);
  if ( greater(projTri, voxelMax-eps) || smaller(projTri-eps, voxelMin) ){
    preciceDebug ( "Found sa in test 2" );
    return;
  }

  // 3. Tests:
  // Triangle edges with voxel sides. 9 (3x3) tests.
  Vector3D edge = vertices[1]; // Edge 0
  edge -= vertices[0];
  Vector3D sepAxis ( 0.0, -edge[2], edge[1] ); // Dim 0: e0 cross edge
  Vector2D projTriVert;
  projTriVert[0] = sepAxis[1]*vertices[0][1] + sepAxis[2]*vertices[0][2];
  projTriVert[1] = sepAxis[1]*vertices[2][1] + sepAxis[2]*vertices[2][2];
  double radius = h[1]*abs(sepAxis[1]) + h[2]*abs(sepAxis[2]);
  if (greater(min(projTriVert), radius-eps) || smaller(max(projTriVert)-eps, -radius)){
    if (greater(radius, 0.0)) return;
  }
  assignList(sepAxis) = edge[2], 0.0, -edge[0]; // Dim 1: e0 cross edge
  projTriVert[0] = sepAxis[0]*vertices[0][0] + sepAxis[2]*vertices[0][2];
  projTriVert[1] = sepAxis[0]*vertices[2][0] + sepAxis[2]*vertices[2][2];
  radius = h[0]*abs(sepAxis[0]) + h[2]*abs(sepAxis[2]);
  if (greater(min(projTriVert), radius-eps) || smaller(max(projTriVert)-eps, -radius)){
    if (greater(radius, 0.0)) return;
  }
  assignList(sepAxis) = -edge[1], edge[0], 0.0; // Dim 2: e0 cross edge
  projTriVert[0] = sepAxis[0]*vertices[0][0] + sepAxis[1]*vertices[0][1];
  projTriVert[1] = sepAxis[0]*vertices[2][0] + sepAxis[1]*vertices[2][1];
  radius = h[0]*abs(sepAxis[0]) + h[1]*abs(sepAxis[1]);
  if (greater(min(projTriVert), radius-eps) || smaller(max(projTriVert)-eps, -radius)){
    if (greater(radius, 0.0)) return;
  }

  edge = vertices[2]; // Edge 1
  edge -= vertices[1];
  assignList(sepAxis) = 0.0, -edge[2], edge[1]; // Dim 0: e1 cross edge
  projTriVert[0] = sepAxis[1]*vertices[0][1] + sepAxis[2]*vertices[0][2];
  projTriVert[1] = sepAxis[1]*vertices[2][1] + sepAxis[2]*vertices[2][2];
  radius = h[1]*abs(sepAxis[1]) + h[2]*abs(sepAxis[2]);
  if (greater(min(projTriVert), radius-eps) || smaller(max(projTriVert)-eps, -radius)){
    if (greater(radius, 0.0)) return;
  }
  assignList(sepAxis) = edge[2], 0.0, -edge[0]; // Dim 1: e1 cross edge
  projTriVert[0] = sepAxis[0]*vertices[0][0] + sepAxis[2]*vertices[0][2];
  projTriVert[1] = sepAxis[0]*vertices[2][0] + sepAxis[2]*vertices[2][2];
  radius = h[0]*abs(sepAxis[0]) + h[2]*abs(sepAxis[2]);
  if (greater(min(projTriVert), radius-eps) || smaller(max(projTriVert)-eps, -radius)){
    if (greater(radius, 0.0)) return;
  }
  assignList(sepAxis) = -edge[1], edge[0], 0.0; // Dim 2: e1 cross edge
  projTriVert[0] = sepAxis[0]*vertices[0][0] + sepAxis[1]*vertices[0][1];
  projTriVert[1] = sepAxis[0]*vertices[2][0] + sepAxis[1]*vertices[2][1];
  radius = h[0]*abs(sepAxis[0]) + h[1]*abs(sepAxis[1]);
  if (greater(min(projTriVert), radius-eps) || smaller(max(projTriVert)-eps, -radius)){
    if (greater(radius, 0.0)) return;
  }

  edge = vertices[0]; // Edge 2
  edge -= vertices[2];
  assignList(sepAxis) = 0.0, -edge[2], edge[1]; // Dim 0: e2 cross edge
  projTriVert[0] = sepAxis[1]*vertices[1][1] + sepAxis[2]*vertices[1][2];
  projTriVert[1] = sepAxis[1]*vertices[2][1] + sepAxis[2]*vertices[2][2];
  radius = h[1]*abs(sepAxis[1]) + h[2]*abs(sepAxis[2]);
  if (greater(min(projTriVert), radius-eps) || smaller(max(projTriVert)-eps, -radius)){
    if (greater(radius, 0.0)) return;
  }
  assignList(sepAxis) = edge[2], 0.0, -edge[0]; // Dim 1: e2 cross edge
  projTriVert[0] = sepAxis[0]*vertices[1][0] + sepAxis[2]*vertices[1][2];
  projTriVert[1] = sepAxis[0]*vertices[2][0] + sepAxis[2]*vertices[2][2];
  radius = h[0]*abs(sepAxis[0]) + h[2]*abs(sepAxis[2]);
  if (greater(min(projTriVert), radius-eps) || smaller(max(projTriVert)-eps, -radius)){
    if (greater(radius, 0.0)) return;
  }
  assignList(sepAxis) = -edge[1], edge[0], 0.0; // Dim 2: e2 cross edge
  projTriVert[0] = sepAxis[0]*vertices[1][0] + sepAxis[1]*vertices[1][1];
  projTriVert[1] = sepAxis[0]*vertices[2][0] + sepAxis[1]*vertices[2][1];
  radius = h[0]*abs(sepAxis[0]) + h[1]*abs(sepAxis[1]);
  if (greater(min(projTriVert), radius-eps) || smaller(max(projTriVert)-eps, -radius)){
    if (greater(radius, 0.0)) return;
  }
//#   ifdef Debug
//    if (_voxelCenter.size() == 3)
//    if (norm2(_voxelCenter - Vector3D(4.343209876543000, 4.0666666666666666, 4.106172839509999)) < 1e-10){
//      precicePrint("Add triangle v0=" << triangle.vertex(0).getCoords()
//                   << ", v1=" << triangle.vertex(1).getCoords()
//                   << ", v2=" << triangle.vertex(2).getCoords());
//    }
//#   endif
  _content.add ( triangle );

# else // not PRECICE_OLD_QUERY --> this section is PRECICE_OLD_QUERY
  typedef utils::GeometryComputations Geocomp;
  bool includeBoundaries = _boundaryInclusion == INCLUDE_BOUNDARY;

  // Test, if circumcircle (circumsphere) of triangle lies completely outside
  Vector3D centerDiff = triangle.getCenter();
  centerDiff -= _voxelCenter;
  abs ( centerDiff, centerDiff );
  if ( oneGreater(centerDiff - triangle.getEnclosingRadius(),  _voxelHalflengths) ){
     return;
  }

   // Test, if triangle center lies completely inside
   Vector3D sidelengths = _voxelHalflengths;
   sidelengths *= 2.0;
   int result = Geocomp::containedInHyperrectangle (
      sidelengths, _voxelCenter, triangle.getCenter() );
   if ( (result == Geocomp::CONTAINED)
        || ((result == Geocomp::TOUCHING) && includeBoundaries) )
   {
      _content.add ( triangle );
      return;
   }

   // Test, if any of the triangle vertices lies inside of the voxel
   result = Geocomp::containedInHyperrectangle (
       sidelengths, _voxelCenter, triangle.vertex(0).getCoords() );
   if ( (result == Geocomp::CONTAINED)
        || ((result == Geocomp::TOUCHING) && includeBoundaries) )
   {
      _content.add ( triangle );
      return;
   }
   result = Geocomp::containedInHyperrectangle (
      sidelengths, _voxelCenter, triangle.vertex(1).getCoords() );
   if ( (result == Geocomp::CONTAINED)
        || ((result == Geocomp::TOUCHING) && includeBoundaries) )
   {
      _content.add ( triangle );
      return;
   }
   result = Geocomp::containedInHyperrectangle (
      sidelengths, _voxelCenter, triangle.vertex(2).getCoords() );
   if ( (result == Geocomp::CONTAINED)
        || ((result == Geocomp::TOUCHING) && includeBoundaries) )
   {
      _content.add ( triangle );
      return;
   }

   // Check for intersections of triangle with voxel edges.
   utils::Vector3D pointOfIntersection;
   utils::Vector3D coords[8];

   // Compute coordinates of corner vertices
   coords[0] = _voxelCenter; coords[0] -= _voxelHalflengths;
   assignList(coords[1]) = _voxelCenter(0) + _voxelHalflengths(0),
                           _voxelCenter(1) - _voxelHalflengths(1),
                           _voxelCenter(2) - _voxelHalflengths(2);
   assignList(coords[2]) = _voxelCenter(0) - _voxelHalflengths(0),
                           _voxelCenter(1) + _voxelHalflengths(1),
                           _voxelCenter(2) - _voxelHalflengths(2);
   assignList(coords[3]) = _voxelCenter(0) + _voxelHalflengths(0),
                           _voxelCenter(1) + _voxelHalflengths(1),
                           _voxelCenter(2) - _voxelHalflengths(2);
   assignList(coords[4]) = _voxelCenter(0) - _voxelHalflengths(0),
                           _voxelCenter(1) - _voxelHalflengths(1),
                           _voxelCenter(2) + _voxelHalflengths(2);
   assignList(coords[5]) = _voxelCenter(0) + _voxelHalflengths(0),
                           _voxelCenter(1) - _voxelHalflengths(1),
                           _voxelCenter(2) + _voxelHalflengths(2);
   assignList(coords[6]) = _voxelCenter(0) - _voxelHalflengths(0),
                           _voxelCenter(1) + _voxelHalflengths(1),
                           _voxelCenter(2) + _voxelHalflengths(2);
   coords[7] = _voxelCenter; coords[7] += _voxelHalflengths;

   // Edges along x-axis

   // Edge 01:
   if ( computeIntersection(triangle, coords[0], coords[1], includeBoundaries) ) {
      _content.add ( triangle );
      return;
   }

   // Edge 23:
   if ( computeIntersection(triangle, coords[2], coords[3], includeBoundaries) ) {
      _content.add ( triangle );
      return;
   }

   // Edge 45:
   if ( computeIntersection(triangle, coords[4], coords[5], includeBoundaries) ) {
      _content.add ( triangle );
      return;
   }

   // Edge 67:
   if ( computeIntersection(triangle, coords[6], coords[7], includeBoundaries) ) {
      _content.add ( triangle );
      return;
   }

   // Edges along y-axis

   // Edge 02:
   if ( computeIntersection(triangle, coords[0], coords[2], includeBoundaries) ) {
      _content.add ( triangle );
      return;
   }

   // Edge 13:
   if ( computeIntersection(triangle, coords[1], coords[3], includeBoundaries) ) {
      _content.add ( triangle );
      return;
   }

   // Edge 46:
   if ( computeIntersection(triangle, coords[4], coords[6], includeBoundaries) ) {
      _content.add ( triangle );
      return;
   }

   // Edge 57:
   if ( computeIntersection(triangle, coords[5], coords[7], includeBoundaries) ) {
      _content.add ( triangle );
      return;
   }

   // Edges along z-axis

   // Edge 04:
   if ( computeIntersection(triangle, coords[0], coords[4], includeBoundaries) ) {
      _content.add ( triangle );
      return;
   }

   // Edge 15:
   if ( computeIntersection(triangle, coords[1], coords[5], includeBoundaries) ) {
      _content.add ( triangle );
      return;
   }

   // Edge 26:
   if ( computeIntersection(triangle, coords[2], coords[6], includeBoundaries) ) {
      _content.add ( triangle );
      return;
   }

   // Edge 37:
   if ( computeIntersection(triangle, coords[3], coords[7], includeBoundaries) ) {
      _content.add ( triangle );
      return;
   }

   // Check if triangle edges interesect voxel sides
   utils::Vector3D center;

   // Side x low
   assignList(center) = _voxelCenter[0] - _voxelHalflengths[0], _voxelCenter[1], _voxelCenter[2];
   coords[0] = triangle.vertex(0).getCoords();
   coords[1] = triangle.vertex(1).getCoords();
   if ( computeIntersection(center, _voxelHalflengths, 0, coords[0], coords[1], includeBoundaries) ) {
     _content.add ( triangle );
     return;
   }
   coords[0] = triangle.vertex(1).getCoords();
   coords[1] = triangle.vertex(2).getCoords();
   if ( computeIntersection(center, _voxelHalflengths, 0, coords[0], coords[1], includeBoundaries) ) {
     _content.add ( triangle );
     return;
   }
   coords[0] = triangle.vertex(2).getCoords();
   coords[1] = triangle.vertex(0).getCoords();
   if ( computeIntersection(center, _voxelHalflengths, 0, coords[0], coords[1], includeBoundaries) ) {
     _content.add ( triangle );
     return;
   }

   // Side x high
   assignList(center) = _voxelCenter[0] + _voxelHalflengths[0], _voxelCenter[1], _voxelCenter[2];
   coords[0] = triangle.vertex(0).getCoords();
   coords[1] = triangle.vertex(1).getCoords();
   if ( computeIntersection(center, _voxelHalflengths, 0, coords[0], coords[1], includeBoundaries) ) {
     _content.add ( triangle );
     return;
   }
   coords[0] = triangle.vertex(1).getCoords();
   coords[1] = triangle.vertex(2).getCoords();
   if ( computeIntersection(center, _voxelHalflengths, 0, coords[0], coords[1], includeBoundaries) ) {
     _content.add ( triangle );
     return;
   }
   coords[0] = triangle.vertex(2).getCoords();
   coords[1] = triangle.vertex(0).getCoords();
   if ( computeIntersection(center, _voxelHalflengths, 0, coords[0], coords[1], includeBoundaries) ) {
     _content.add ( triangle );
     return;
   }

   // Side y low
   assignList(center) = _voxelCenter[0], _voxelCenter[1] - _voxelHalflengths[1], _voxelCenter[2];
   coords[0] = triangle.vertex(0).getCoords();
   coords[1] = triangle.vertex(1).getCoords();
   if ( computeIntersection(center, _voxelHalflengths, 1, coords[0], coords[1], includeBoundaries) ) {
     _content.add ( triangle );
     return;
   }
   coords[0] = triangle.vertex(1).getCoords();
   coords[1] = triangle.vertex(2).getCoords();
   if ( computeIntersection(center, _voxelHalflengths, 1, coords[0], coords[1], includeBoundaries) ) {
     _content.add ( triangle );
     return;
   }
   coords[0] = triangle.vertex(2).getCoords();
   coords[1] = triangle.vertex(0).getCoords();
   if ( computeIntersection(center, _voxelHalflengths, 1, coords[0], coords[1], includeBoundaries) ) {
     _content.add ( triangle );
     return;
   }

   // Side y high
   assignList(center) = _voxelCenter[0], _voxelCenter[1] + _voxelHalflengths[1], _voxelCenter[2];
   coords[0] = triangle.vertex(0).getCoords();
   coords[1] = triangle.vertex(1).getCoords();
   if ( computeIntersection(center, _voxelHalflengths, 1, coords[0], coords[1], includeBoundaries) ) {
     _content.add ( triangle );
     return;
   }
   coords[0] = triangle.vertex(1).getCoords();
   coords[1] = triangle.vertex(2).getCoords();
   if ( computeIntersection(center, _voxelHalflengths, 1, coords[0], coords[1], includeBoundaries) ) {
     _content.add ( triangle );
     return;
   }
   coords[0] = triangle.vertex(2).getCoords();
   coords[1] = triangle.vertex(0).getCoords();
   if ( computeIntersection(center, _voxelHalflengths, 1, coords[0], coords[1], includeBoundaries) ) {
     _content.add ( triangle );
     return;
   }

   // Side z low
   assignList(center) = _voxelCenter[0], _voxelCenter[1], _voxelCenter[2] - _voxelHalflengths[2];
   coords[0] = triangle.vertex(0).getCoords();
   coords[1] = triangle.vertex(1).getCoords();
   if ( computeIntersection(center, _voxelHalflengths, 2, coords[0], coords[1], includeBoundaries) ) {
     _content.add ( triangle );
     return;
   }
   coords[0] = triangle.vertex(1).getCoords();
   coords[1] = triangle.vertex(2).getCoords();
   if ( computeIntersection(center, _voxelHalflengths, 2, coords[0], coords[1], includeBoundaries) ) {
     _content.add ( triangle );
     return;
   }
   coords[0] = triangle.vertex(2).getCoords();
   coords[1] = triangle.vertex(0).getCoords();
   if ( computeIntersection(center, _voxelHalflengths, 2, coords[0], coords[1], includeBoundaries) ) {
     _content.add ( triangle );
     return;
   }

   // Side z high
   assignList(center) = _voxelCenter[0], _voxelCenter[1], _voxelCenter[2] + _voxelHalflengths[2];
   coords[0] = triangle.vertex(0).getCoords();
   coords[1] = triangle.vertex(1).getCoords();
   if ( computeIntersection(center, _voxelHalflengths, 2, coords[0], coords[1], includeBoundaries) ) {
     _content.add ( triangle );
     return;
   }
   coords[0] = triangle.vertex(1).getCoords();
   coords[1] = triangle.vertex(2).getCoords();
   if ( computeIntersection(center, _voxelHalflengths, 2, coords[0], coords[1], includeBoundaries) ) {
     _content.add ( triangle );
     return;
   }
   coords[0] = triangle.vertex(2).getCoords();
   coords[1] = triangle.vertex(0).getCoords();
   if ( computeIntersection(center, _voxelHalflengths, 2, coords[0], coords[1], includeBoundaries) ) {
     _content.add ( triangle );
     return;
   }

   // Check touching - intersecting triangles
   //std::vector<Vector3D> coords;
   coords[0] = triangle.vertex(0).getCoords();
   coords[1] = triangle.vertex(1).getCoords();
   coords[2] = triangle.vertex(2).getCoords();
   coords[3] = triangle.vertex(0).getCoords();
   for ( int iVertex=0; iVertex < 3; iVertex++ ) {
     Vector3D vertex0 ( coords[iVertex] );
     Vector3D vertex1 ( coords[iVertex+1] );
     Vector3D toVertex ( _voxelCenter );
     toVertex -= vertex0;
     abs ( toVertex, toVertex );
     bool touching = not oneGreater ( toVertex, _voxelHalflengths );
     touching &= not allGreater ( _voxelHalflengths, toVertex );
     typedef utils::GeometryComputations GeoComp;
     Vector3D edgeDirection ( vertex1 );
     edgeDirection -= vertex0;
     edgeDirection /= norm2 ( edgeDirection );
     edgeDirection *= 10.0 * NUMERICAL_ZERO_DIFFERENCE;
     toVertex = _voxelCenter;
     toVertex -= vertex0;
     toVertex -= edgeDirection;
     abs ( toVertex, toVertex );
     if ( allGreater(_voxelHalflengths, toVertex) ) {
       _content.add ( triangle );
       return;
     }
   }

   for ( int iVertex=0; iVertex < 3; iVertex++ ) {
     Vector3D vertex ( triangle.vertex(iVertex).getCoords() );
     Vector3D toVertex ( _voxelCenter );
     toVertex -= vertex;
     abs ( toVertex, toVertex );
     bool touching = not oneGreater ( toVertex, _voxelHalflengths );
     touching &= not allGreater ( _voxelHalflengths, toVertex );
     Vector3D direction ( triangle.getCenter() );
     direction -= vertex;
     direction /= norm2 ( direction );
     direction *= 10.0 * NUMERICAL_ZERO_DIFFERENCE;
     toVertex = _voxelCenter;
     toVertex -= vertex;
     toVertex -= direction;
     abs ( toVertex, toVertex );
     if ( allGreater(_voxelHalflengths, toVertex) ) {
       _content.add ( triangle );
       return;
     }
   }
# endif // PRECICE_OLD_QUERY
}

bool FindVoxelContent:: computeIntersection
(
  const utils::Vector3D&  squareCenter,
  const utils::DynVector& halflengths,
  int                     squareNormalDirection,
  const utils::Vector3D&  firstPointSegment,
  const utils::Vector3D&  secondPointSegment,
  bool                    countTouchingAsIntersection ) const
{
  preciceTrace ( "computeIntersection()", squareCenter, halflengths,
                  squareNormalDirection, firstPointSegment, secondPointSegment,
                  countTouchingAsIntersection );
  using namespace tarch::la;
  assertion ( (squareNormalDirection >= 0) && (squareNormalDirection < 3) );
  typedef utils::GeometryComputations GeoComp;
  using utils::Vector2D;
  using utils::Vector3D;
  Vector3D normal ( 0.0 );
  Vector3D intersection ( 0.0 );
  normal[squareNormalDirection] = 1.0;
  int result = GeoComp::segmentPlaneIntersection (
      squareCenter, normal, firstPointSegment, secondPointSegment, intersection );
  if ( result == GeoComp::NO_INTERSECTION ) {
//    precicePrint ( "computeIntersection(): no square plane intersection" );
    return false;
  }
  else if ( result == GeoComp::CONTAINED ) {
    if ( not countTouchingAsIntersection ) {
      return false;
    }
    // Edge is contained in plane of square.
    // Do first check, whether one of the points of the edge is contained in
    // the square. If not, check, if any of the squares's edges does
    // intersect with the segment.
    int indices[2];
    if (squareNormalDirection == 0){
      indices[0] = 1;
      indices[1] = 2;
    }
    else if (squareNormalDirection == 1){
      indices[0] = 0;
      indices[1] = 2;
    }
    else {
      assertion ( squareNormalDirection == 2, squareNormalDirection );
      indices[0] = 0;
      indices[1] = 1;
    }
    Vector2D center2D ( squareCenter[indices[0]], squareCenter[indices[1]] );
    Vector2D q ( firstPointSegment[indices[0]], firstPointSegment[indices[1]] );
    Vector2D r ( secondPointSegment[indices[0]], secondPointSegment[indices[1]] );
    Vector2D halflengths2D ( halflengths[indices[0]], halflengths[indices[1]] );
    Vector2D qcenter ( q - center2D );
    abs ( qcenter, qcenter );
    Vector2D rcenter ( r - center2D );
    abs ( rcenter, rcenter );
    bool qOutside = oneGreater ( qcenter, halflengths2D );
    bool qInside = allGreater ( halflengths2D, qcenter );
    bool rOutside = oneGreater ( rcenter, halflengths2D );
    bool rInside = allGreater ( halflengths2D, rcenter );

    if ( qInside || rInside ) { // One or both segment points lie inside
//      precicePrint ( "q and/or r lie in square (segment in plane of square)" );
      return true;
    }
    else if ( (not (qOutside && rOutside)) && countTouchingAsIntersection ) {
      // One or both segment points touch the square
//      precicePrint ( "q and/or r touch square (segment in plane of square)" );
      return true;
    }
    else { // No point is inside, the segment might interesect still
      Vector2D a ( center2D[0] - halflengths2D[0], center2D[1] - halflengths2D[1] );
      Vector2D b ( center2D[0] + halflengths2D[0], center2D[1] - halflengths2D[1] );
      Vector2D c ( center2D[0] + halflengths2D[0], center2D[1] + halflengths2D[1] );
      Vector2D d ( center2D[0] - halflengths2D[0], center2D[1] + halflengths2D[1] );
      bool intersect = false;
      intersect |= GeoComp::segmentsIntersect(a, b, q, r, countTouchingAsIntersection);
      intersect |= GeoComp::segmentsIntersect(b, c, q, r, countTouchingAsIntersection);
      intersect |= GeoComp::segmentsIntersect(c, d, q, r, countTouchingAsIntersection);
      intersect |= GeoComp::segmentsIntersect(d, a, q, r, countTouchingAsIntersection);
      return intersect;
    }
  }
  else if ( (result == GeoComp::TOUCHING) && (! countTouchingAsIntersection) ) {
    return false;
  }
  assertion ( (result == GeoComp::INTERSECTION) || (result == GeoComp::TOUCHING));
  // Segment intersects or touches plane of square, see if intersection point is
  // contained or touches square
  int indices[2];
  if (squareNormalDirection == 0){
    indices[0] = 1;
    indices[1] = 2;
  }
  else if (squareNormalDirection == 1){
    indices[0] = 0;
    indices[1] = 2;
  }
  else {
    assertion ( squareNormalDirection == 3, squareNormalDirection );
    indices[0] = 0;
    indices[1] = 1;
  }
  Vector2D center2D ( squareCenter[indices[0]], squareCenter[indices[1]] );
  Vector2D intersection2D ( intersection[indices[0]], intersection[indices[1]] );
  Vector2D halflengths2D ( halflengths[indices[0]], halflengths[indices[1]] );
  Vector2D toCenter ( intersection2D - center2D );
  abs ( toCenter, toCenter );
  bool outside = oneGreater ( toCenter, halflengths2D );
  bool inside = allGreater ( halflengths2D, toCenter );

  if ( inside ) {
    return true;
  }
  else if ( (not outside) && countTouchingAsIntersection  ) {
    return true;
  }
  return false;
}

bool FindVoxelContent:: computeIntersection
(
  const mesh::Triangle&  triangle,
  const utils::Vector3D& firstPointEdge,
  const utils::Vector3D& secondPointEdge,
  bool                   countTouchingAsIntersection )
{
  preciceTrace ( "computeIntersection()", triangle.vertex(0).getCoords(),
                  triangle.vertex(1).getCoords(), triangle.vertex(2).getCoords(),
                  firstPointEdge, secondPointEdge, countTouchingAsIntersection );
  typedef utils::GeometryComputations Geocomp;
  using utils::Vector2D;
  using utils::Vector3D;
  using namespace tarch::la;
  Vector3D pointOfIntersection;
  int result = Geocomp::segmentPlaneIntersection (
    triangle.edge(0).getCenter(), triangle.getNormal(),
    firstPointEdge, secondPointEdge, pointOfIntersection );
  if ( result == Geocomp::NO_INTERSECTION ) {
    return false;
  }
  else if ( result == Geocomp::CONTAINED ) {
    if ( not countTouchingAsIntersection ) {
      return false;
    }
    // Edge is contained in plane of triangle.
    // Do first check, whether one of the points of the edge is contained in
    // the triangle. If not, check, if any of the triangle's edges does
    // intersect with the segment.
    Vector3D normal ( triangle.getNormal() );
    int indexMax = tarch::la::indexMax ( abs(normal, normal) );
    Vector2D a (Geocomp::projectVector(triangle.vertex(0).getCoords(), indexMax));
    Vector2D b (Geocomp::projectVector(triangle.vertex(1).getCoords(), indexMax));
    Vector2D c (Geocomp::projectVector(triangle.vertex(2).getCoords(), indexMax));
    Vector2D q (Geocomp::projectVector(firstPointEdge, indexMax));
    int containedResult = Geocomp::containedInTriangle(a, b, c, q);

    if  ( (containedResult == Geocomp::TOUCHING) &&  countTouchingAsIntersection ){
      return true;
    }
    Vector2D r (Geocomp::projectVector(secondPointEdge, indexMax));
    containedResult = Geocomp::containedInTriangle(a, b, c, r);
    if  ( (containedResult == Geocomp::TOUCHING) && countTouchingAsIntersection ){
      return true;
    }
    if ( Geocomp::segmentsIntersect(a, b, q, r, countTouchingAsIntersection) ||
         Geocomp::segmentsIntersect(b, c, q, r, countTouchingAsIntersection) ||
         Geocomp::segmentsIntersect(c, a, q, r, countTouchingAsIntersection) )
    {
      return true;
    }
    return false;
  }
  else if ( (result == Geocomp::TOUCHING) && (not countTouchingAsIntersection) ) {
    return false;
  }

  assertion ((result == Geocomp::INTERSECTION) || (result == Geocomp::TOUCHING));
  // Compute signed areas of triangle and intersection point
  Vector3D normal ( triangle.getNormal() );
  int indexMin = tarch::la::indexMax ( abs(normal,normal) );
  result = Geocomp::containedInTriangle (
      Geocomp::projectVector(triangle.vertex(0).getCoords(), indexMin),
      Geocomp::projectVector(triangle.vertex(1).getCoords(), indexMin),
      Geocomp::projectVector(triangle.vertex(2).getCoords(), indexMin),
      Geocomp::projectVector ( pointOfIntersection, indexMin ) );
  if ( result == Geocomp::NOT_CONTAINED ) {
    return false;
  }
  else if ( (result == Geocomp::TOUCHING) && (not countTouchingAsIntersection) ) {
    return false;
  }
  return true;
}

}} // namespace precice, query
