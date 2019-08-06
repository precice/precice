#include "FindVoxelContent.hpp"
#include "math/geometry.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Mesh.hpp"
#include "math/math.hpp"

namespace precice {
namespace query {

FindVoxelContent:: FindVoxelContent
(
  const Eigen::VectorXd&  voxelCenter,
  const Eigen::VectorXd&  halflengths,
  BoundaryInclusion boundaryInclusion )
:
  _voxelCenter ( voxelCenter ),
  _voxelHalflengths ( halflengths ),
  _boundaryInclusion ( boundaryInclusion ),
  _dimensions ( voxelCenter.size() )
{
  P_TRACE(voxelCenter, halflengths, boundaryInclusion);
  P_ASSERT( voxelCenter.size() == halflengths.size(),
              voxelCenter.size(), halflengths.size() );
  P_ASSERT( (_dimensions == 2) || (_dimensions == 3), _dimensions );
}

const Eigen::VectorXd& FindVoxelContent:: getVoxelCenter() const
{
  return _voxelCenter;
}

const Eigen::VectorXd& FindVoxelContent:: getVoxelHalflengths() const
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

void FindVoxelContent::checkVertex( mesh::Vertex& vertex )
{
  P_TRACE();
  P_ASSERT(vertex.getDimensions() == _dimensions, vertex.getDimensions(),
            _dimensions);
  Eigen::VectorXd toVertex = vertex.getCoords();
  toVertex -= _voxelCenter;
  toVertex = toVertex.cwiseAbs();
  Eigen::VectorXd absHalflengths = Eigen::VectorXd::Zero(_voxelHalflengths.size());
  absHalflengths = _voxelHalflengths.cwiseAbs();
  if (_boundaryInclusion == INCLUDE_BOUNDARY) {
    if (math::oneGreater(toVertex, absHalflengths)) 
      return;
  }
  else
    if (math::oneGreaterEquals(toVertex, absHalflengths))
      return;
  
  _content.add(vertex);
}

void FindVoxelContent::checkEdge( mesh::Edge& edge )
{
  P_TRACE(edge.vertex(0).getCoords(), edge.vertex(1).getCoords());
  using Eigen::Vector3d; using Eigen::Vector2d; using Eigen::VectorXd;
  using math::greater; using math::smaller;
  P_ASSERT(edge.getDimensions() == _dimensions, edge.getDimensions(),
            _dimensions);

  // For excluding boundary case, add eps to achieve greater/smaller-equals
  double eps = 0.0;
  if ( _boundaryInclusion == EXCLUDE_BOUNDARY ){
    eps = 2.0 * math::NUMERICAL_ZERO_DIFFERENCE;
  }

  if ( _dimensions == 2 ){
    // Test if edge intersects with rectangle, using seperating axis theorem
    Vector2d a = edge.vertex(0).getCoords();
    a -= _voxelCenter; // Move to origin
    Vector2d b = edge.vertex(1).getCoords();
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
    VectorXd edgeNormal = edge.getNormal();
    int dim;
    edgeNormal.cwiseAbs().maxCoeff(&dim);  // 1D dimensions

    double projEdge = edgeNormal.dot(a) * edgeNormal[dim];
    Eigen::Vector4d projVoxel;
    Vector2d corner ( -_voxelHalflengths[0], -_voxelHalflengths[1] );
    projVoxel[0] = edgeNormal.dot(corner) * edgeNormal[dim];
    corner << _voxelHalflengths[0], -_voxelHalflengths[1];
    projVoxel[1] = edgeNormal.dot(corner) * edgeNormal[dim];
    corner << -_voxelHalflengths[0], _voxelHalflengths[1];
    projVoxel[2] = edgeNormal.dot(corner) * edgeNormal[dim];
    corner << _voxelHalflengths[0], _voxelHalflengths[1];
    projVoxel[3] = edgeNormal.dot(corner) * edgeNormal[dim];
    double voxelMax = projVoxel.maxCoeff();
    double voxelMin = projVoxel.minCoeff();
    if ( greater(projEdge, voxelMax-eps) || smaller(projEdge-eps, voxelMin) ){
      return;
    }
    // If all tests pass, the edge intersects the voxel
    _content.add (edge);

  }
  else { // 3D
    P_ASSERT(_dimensions == 3, _dimensions);
    // Test if edge intersects with rectangle, using seperating axis theorem
    Vector3d a = edge.vertex(0).getCoords();
    a -= _voxelCenter; // Move to origin
    Vector3d b = edge.vertex(1).getCoords();
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

    Vector3d edgeVector = b - a;
    // Separation axis is formed by cross product of unit vector in x-direction
    // and edge vector:
    Vector3d sepAxis(0.0, -edgeVector[2], edgeVector[1]); // Dim 0
    // Edge is projected onto separation axis by dot product:
    double proj = sepAxis[1]*a[1] + sepAxis[2]*a[2];
    double radius = _voxelHalflengths[1]*std::abs(sepAxis[1]) + _voxelHalflengths[2]*std::abs(sepAxis[2]);
    if (greater(proj, radius-eps) || smaller(proj-eps, -radius)){
      if (radius > 0.0) return;
    }
    sepAxis << edgeVector[2], 0.0, -edgeVector[0]; // Dim 1
    proj = sepAxis[0]*a[0] + sepAxis[2]*a[2];
    radius = _voxelHalflengths[0]*std::abs(sepAxis[0]) + _voxelHalflengths[2]*std::abs(sepAxis[2]);
    if (greater(proj, radius-eps) || smaller(proj-eps, -radius)){
      if (radius > 0.0) return;
    }
    sepAxis << -edgeVector[1], edgeVector[0], 0.0; // Dim 2
    proj = sepAxis[0]*a[0] + sepAxis[1]*a[1];
    radius = _voxelHalflengths[0]*std::abs(sepAxis[0]) + _voxelHalflengths[1]*std::abs(sepAxis[1]);
    if (greater(proj, radius-eps) || smaller(proj-eps, -radius)){
      if (radius > 0.0) return;
    }

    // If all tests pass, the edge intersects the voxel
    _content.add(edge);
  }
}

void FindVoxelContent:: checkTriangle
(
  mesh::Triangle& triangle )
{
  P_TRACE(triangle.getID(), triangle.getCenter());
  P_ASSERT( _dimensions == 3, _dimensions );
  using math::greater; using math::smaller;
  using Eigen::Vector3d;
  using Eigen::Vector2d;
  using std::abs;

  // Use Separating axis theorem to find intersection

  // For excluding boundary case, add eps to achieve greater/smaller-equals
  double eps = 0.0;
  if ( _boundaryInclusion == EXCLUDE_BOUNDARY ){
    eps = 2.0 * math::NUMERICAL_ZERO_DIFFERENCE;
  }

  // Shift triangle vertex coordinates to have voxel at coordinate origin
  std::vector<Vector3d> vertices(3);
  for (int i=0; i < 3; i++){
    vertices[i] = triangle.vertex(i).getCoords();
    vertices[i] -= _voxelCenter;
  }

  // 1. Tests:
  // Voxel sides against minimal triangle axis-aligned bounding box
  Vector3d coords; // Holds triangle coords of one dimension
  double triangleMin = 0.0;
  double triangleMax = 0.0;
  double voxelMax = 0.0;
  double voxelMin = 0.0;
  for (int dim=0; dim < 3; dim++ ){     // Test projection in "dim"
    for (int i=0; i < 3; i++){
      coords[i] = vertices[i][dim];
    }
    triangleMin = coords.minCoeff();
    voxelMax = _voxelHalflengths[dim];
    if ( greater(triangleMin, voxelMax-eps) ){
      P_DEBUG( "Found sa in test 1.1" );
      return;
    }
    triangleMax = coords.maxCoeff();
    voxelMin = -1.0 * _voxelHalflengths[dim];
    if ( smaller(triangleMax-eps, voxelMin) ){
      P_DEBUG( "Found sa in test 1.2" );
      return;
    }
  }

  // 2. Tests:
  // Triangle plane intersecting with voxel. Projection of voxel vertices on
  // triangle normal.
  Vector3d n = triangle.getNormal();
  int dim = -1; // 1D dimensions
  n.cwiseAbs().maxCoeff(&dim);  // 1D dimensions

  double projTri = n.dot(vertices[0]) * n[dim];
  Eigen::Matrix<double, 8, 1> projVoxel;
  Eigen::VectorXd& h = _voxelHalflengths;
  Vector3d corner ( -h[0], -h[1], -h[2] );
  projVoxel[0] = n.dot(corner) * n[dim];
  corner <<  h[0], -h[1], -h[2];
  projVoxel[1] = n.dot(corner) * n[dim];
  corner << -h[0],  h[1], -h[2];
  projVoxel[2] = n.dot(corner) * n[dim];
  corner <<  h[0],  h[1], -h[2];
  projVoxel[3] = n.dot(corner) * n[dim];
  corner << -h[0], -h[1],  h[2];
  projVoxel[4] = n.dot(corner) * n[dim];
  corner <<  h[0], -h[1],  h[2];
  projVoxel[5] = n.dot(corner) * n[dim];
  corner << -h[0],  h[1],  h[2];
  projVoxel[6] = n.dot(corner) * n[dim];
  corner <<  h[0],  h[1],  h[2];
  projVoxel[7] = n.dot(corner) * n[dim];
  voxelMax = projVoxel.maxCoeff();
  voxelMin = projVoxel.minCoeff();
  if ( greater(projTri, voxelMax-eps) || smaller(projTri-eps, voxelMin) ){
    P_DEBUG( "Found sa in test 2" );
    return;
  }

  // 3. Tests:
  // Triangle edges with voxel sides. 9 (3x3) tests.
  Vector3d edge = vertices[1]; // Edge 0
  edge -= vertices[0];
  Vector3d sepAxis ( 0.0, -edge[2], edge[1] ); // Dim 0: e0 cross edge
  Vector2d projTriVert;
  projTriVert[0] = sepAxis[1]*vertices[0][1] + sepAxis[2]*vertices[0][2];
  projTriVert[1] = sepAxis[1]*vertices[2][1] + sepAxis[2]*vertices[2][2];
  double radius = h[1]*abs(sepAxis[1]) + h[2]*abs(sepAxis[2]);
  if (greater(projTriVert.minCoeff(), radius-eps) || smaller(projTriVert.maxCoeff()-eps, -radius)){
    if (greater(radius, 0.0)) return;
  }
  sepAxis << edge[2], 0.0, -edge[0]; // Dim 1: e0 cross edge
  projTriVert[0] = sepAxis[0]*vertices[0][0] + sepAxis[2]*vertices[0][2];
  projTriVert[1] = sepAxis[0]*vertices[2][0] + sepAxis[2]*vertices[2][2];
  radius = h[0]*abs(sepAxis[0]) + h[2]*abs(sepAxis[2]);
  if (greater(projTriVert.minCoeff(), radius-eps) || smaller(projTriVert.maxCoeff()-eps, -radius)){
    if (greater(radius, 0.0)) return;
  }
  sepAxis << -edge[1], edge[0], 0.0; // Dim 2: e0 cross edge
  projTriVert[0] = sepAxis[0]*vertices[0][0] + sepAxis[1]*vertices[0][1];
  projTriVert[1] = sepAxis[0]*vertices[2][0] + sepAxis[1]*vertices[2][1];
  radius = h[0]*abs(sepAxis[0]) + h[1]*abs(sepAxis[1]);
  if (greater(projTriVert.minCoeff(), radius-eps) || smaller(projTriVert.maxCoeff()-eps, -radius)) {
    if (greater(radius, 0.0)) return;
  }

  edge = vertices[2]; // Edge 1
  edge -= vertices[1];
  sepAxis << 0.0, -edge[2], edge[1]; // Dim 0: e1 cross edge
  projTriVert[0] = sepAxis[1]*vertices[0][1] + sepAxis[2]*vertices[0][2];
  projTriVert[1] = sepAxis[1]*vertices[2][1] + sepAxis[2]*vertices[2][2];
  radius = h[1]*abs(sepAxis[1]) + h[2]*abs(sepAxis[2]);
  if (greater(projTriVert.minCoeff(), radius-eps) || smaller(projTriVert.maxCoeff()-eps, -radius)){
    if (greater(radius, 0.0)) return;
  }
  sepAxis << edge[2], 0.0, -edge[0]; // Dim 1: e1 cross edge
  projTriVert[0] = sepAxis[0]*vertices[0][0] + sepAxis[2]*vertices[0][2];
  projTriVert[1] = sepAxis[0]*vertices[2][0] + sepAxis[2]*vertices[2][2];
  radius = h[0]*abs(sepAxis[0]) + h[2]*abs(sepAxis[2]);
  if (greater(projTriVert.minCoeff(), radius-eps) || smaller(projTriVert.maxCoeff()-eps, -radius)) {
    if (greater(radius, 0.0)) return;
  }
  sepAxis << -edge[1], edge[0], 0.0; // Dim 2: e1 cross edge
  projTriVert[0] = sepAxis[0]*vertices[0][0] + sepAxis[1]*vertices[0][1];
  projTriVert[1] = sepAxis[0]*vertices[2][0] + sepAxis[1]*vertices[2][1];
  radius = h[0]*abs(sepAxis[0]) + h[1]*abs(sepAxis[1]);
  if (greater(projTriVert.minCoeff(), radius-eps) || smaller(projTriVert.maxCoeff()-eps, -radius)){
    if (greater(radius, 0.0)) return;
  }

  edge = vertices[0]; // Edge 2
  edge -= vertices[2];
  sepAxis << 0.0, -edge[2], edge[1]; // Dim 0: e2 cross edge
  projTriVert[0] = sepAxis[1]*vertices[1][1] + sepAxis[2]*vertices[1][2];
  projTriVert[1] = sepAxis[1]*vertices[2][1] + sepAxis[2]*vertices[2][2];
  radius = h[1]*abs(sepAxis[1]) + h[2]*abs(sepAxis[2]);
  if (greater(projTriVert.minCoeff(), radius-eps) || smaller(projTriVert.maxCoeff()-eps, -radius)) {
    if (greater(radius, 0.0)) return;
  }
  sepAxis << edge[2], 0.0, -edge[0]; // Dim 1: e2 cross edge
  projTriVert[0] = sepAxis[0]*vertices[1][0] + sepAxis[2]*vertices[1][2];
  projTriVert[1] = sepAxis[0]*vertices[2][0] + sepAxis[2]*vertices[2][2];
  radius = h[0]*abs(sepAxis[0]) + h[2]*abs(sepAxis[2]);
  if (greater(projTriVert.minCoeff(), radius-eps) || smaller(projTriVert.maxCoeff()-eps, -radius)){
    if (greater(radius, 0.0)) return;
  }
  sepAxis << -edge[1], edge[0], 0.0; // Dim 2: e2 cross edge
  projTriVert[0] = sepAxis[0]*vertices[1][0] + sepAxis[1]*vertices[1][1];
  projTriVert[1] = sepAxis[0]*vertices[2][0] + sepAxis[1]*vertices[2][1];
  radius = h[0]*abs(sepAxis[0]) + h[1]*abs(sepAxis[1]);
  if (greater(projTriVert.minCoeff(), radius-eps) || smaller(projTriVert.maxCoeff()-eps, -radius)){
    if (greater(radius, 0.0)) return;
  }
  _content.add ( triangle );
}

bool FindVoxelContent:: computeIntersection
(
  const Eigen::Vector3d& squareCenter,
  const Eigen::VectorXd& halflengths,
  int                    squareNormalDirection,
  const Eigen::Vector3d& firstPointSegment,
  const Eigen::Vector3d& secondPointSegment,
  bool                   countTouchingAsIntersection ) const
{
  P_TRACE(squareCenter, halflengths, squareNormalDirection,
        firstPointSegment, secondPointSegment, countTouchingAsIntersection);
  
  P_ASSERT( (squareNormalDirection >= 0) && (squareNormalDirection < 3) );
  namespace geo = math::geometry;
  using Eigen::Vector2d;
  using Eigen::Vector3d;
  Vector3d normal = Vector3d::Zero();
  Vector3d intersection = Vector3d::Zero();
  normal[squareNormalDirection] = 1.0;
  int result = geo::segmentPlaneIntersection (
      squareCenter, normal, firstPointSegment, secondPointSegment, intersection );
  if ( result == geo::NO_INTERSECTION ) {
//    P_INFO( "computeIntersection(): no square plane intersection" );
    return false;
  }
  else if ( result == geo::CONTAINED ) {
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
      P_ASSERT( squareNormalDirection == 2, squareNormalDirection );
      indices[0] = 0;
      indices[1] = 1;
    }
    Vector2d center2D ( squareCenter[indices[0]], squareCenter[indices[1]] );
    Vector2d q ( firstPointSegment[indices[0]], firstPointSegment[indices[1]] );
    Vector2d r ( secondPointSegment[indices[0]], secondPointSegment[indices[1]] );
    Vector2d halflengths2D ( halflengths[indices[0]], halflengths[indices[1]] );
    Vector2d qcenter = q - center2D;
    qcenter = qcenter.cwiseAbs();
    Vector2d rcenter = r - center2D;
    rcenter = rcenter.cwiseAbs();
    bool qOutside = math::oneGreater ( qcenter, halflengths2D );
    bool qInside = math::allGreater ( halflengths2D, qcenter );
    bool rOutside = math::oneGreater ( rcenter, halflengths2D );
    bool rInside = math::allGreater ( halflengths2D, rcenter );
    
    if ( qInside || rInside ) { // One or both segment points lie inside
//      P_INFO( "q and/or r lie in square (segment in plane of square)" );
      return true;
    }
    else if ( (not (qOutside && rOutside)) && countTouchingAsIntersection ) {
      // One or both segment points touch the square
//      P_INFO( "q and/or r touch square (segment in plane of square)" );
      return true;
    }
    else { // No point is inside, the segment might interesect still
      Vector2d a ( center2D[0] - halflengths2D[0], center2D[1] - halflengths2D[1] );
      Vector2d b ( center2D[0] + halflengths2D[0], center2D[1] - halflengths2D[1] );
      Vector2d c ( center2D[0] + halflengths2D[0], center2D[1] + halflengths2D[1] );
      Vector2d d ( center2D[0] - halflengths2D[0], center2D[1] + halflengths2D[1] );
      bool intersect = false;
      intersect |= geo::segmentsIntersect(a, b, q, r, countTouchingAsIntersection);
      intersect |= geo::segmentsIntersect(b, c, q, r, countTouchingAsIntersection);
      intersect |= geo::segmentsIntersect(c, d, q, r, countTouchingAsIntersection);
      intersect |= geo::segmentsIntersect(d, a, q, r, countTouchingAsIntersection);
      return intersect;
    }
  }
  else if ( (result == geo::TOUCHING) && (! countTouchingAsIntersection) ) {
    return false;
  }
  P_ASSERT( (result == geo::INTERSECTION) || (result == geo::TOUCHING));
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
    P_ASSERT( squareNormalDirection == 3, squareNormalDirection );
    indices[0] = 0;
    indices[1] = 1;
  }
  Vector2d center2D ( squareCenter[indices[0]], squareCenter[indices[1]] );
  Vector2d intersection2D ( intersection[indices[0]], intersection[indices[1]] );
  Vector2d halflengths2D ( halflengths[indices[0]], halflengths[indices[1]] );
  Vector2d toCenter = intersection2D - center2D;
  toCenter = toCenter.cwiseAbs();
  bool outside = math::oneGreater ( toCenter, halflengths2D );
  bool inside = math::allGreater ( halflengths2D, toCenter );
  
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
  const Eigen::Vector3d& firstPointEdge,
  const Eigen::Vector3d& secondPointEdge,
  bool                   countTouchingAsIntersection )
{
  P_TRACE(triangle.vertex(0).getCoords(), triangle.vertex(1).getCoords(), triangle.vertex(2).getCoords(),
        firstPointEdge, secondPointEdge, countTouchingAsIntersection );
  namespace geo = math::geometry;
  using Eigen::Vector2d;
  using Eigen::Vector3d;
  Vector3d pointOfIntersection;
  int result = geo::segmentPlaneIntersection (
    triangle.edge(0).getCenter(), triangle.getNormal(),
    firstPointEdge, secondPointEdge, pointOfIntersection );
  if ( result == geo::NO_INTERSECTION ) {
    return false;
  }
  else if ( result == geo::CONTAINED ) {
    if ( not countTouchingAsIntersection ) {
      return false;
    }
    // Edge is contained in plane of triangle.
    // Do first check, whether one of the points of the edge is contained in
    // the triangle. If not, check, if any of the triangle's edges does
    // intersect with the segment.
    Vector3d normal ( triangle.getNormal() );
    int indexMax;
    normal.cwiseAbs().maxCoeff(&indexMax);
    Vector2d a (geo::projectVector(triangle.vertex(0).getCoords(), indexMax));
    Vector2d b (geo::projectVector(triangle.vertex(1).getCoords(), indexMax));
    Vector2d c (geo::projectVector(triangle.vertex(2).getCoords(), indexMax));
    Vector2d q (geo::projectVector(firstPointEdge, indexMax));
    int containedResult = geo::containedInTriangle(a, b, c, q);

    if  ( (containedResult == geo::TOUCHING) &&  countTouchingAsIntersection ){
      return true;
    }
    Vector2d r (geo::projectVector(secondPointEdge, indexMax));
    containedResult = geo::containedInTriangle(a, b, c, r);
    if  ( (containedResult == geo::TOUCHING) && countTouchingAsIntersection ){
      return true;
    }
    if ( geo::segmentsIntersect(a, b, q, r, countTouchingAsIntersection) ||
         geo::segmentsIntersect(b, c, q, r, countTouchingAsIntersection) ||
         geo::segmentsIntersect(c, a, q, r, countTouchingAsIntersection) )
    {
      return true;
    }
    return false;
  }
  else if ( (result == geo::TOUCHING) && (not countTouchingAsIntersection) ) {
    return false;
  }

  P_ASSERT((result == geo::INTERSECTION) || (result == geo::TOUCHING));
  // Compute signed areas of triangle and intersection point
  Vector3d normal ( triangle.getNormal() );
  int indexMin;
  normal.cwiseAbs().maxCoeff(&indexMin);
  result = geo::containedInTriangle (
      geo::projectVector(triangle.vertex(0).getCoords(), indexMin),
      geo::projectVector(triangle.vertex(1).getCoords(), indexMin),
      geo::projectVector(triangle.vertex(2).getCoords(), indexMin),
      geo::projectVector ( pointOfIntersection, indexMin ) );
  if ( result == geo::NOT_CONTAINED ) {
    return false;
  }
  else if ( (result == geo::TOUCHING) && (not countTouchingAsIntersection) ) {
    return false;
  }
  return true;
}

}} // namespace precice, query
