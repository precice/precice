#include "BalanceVertexPositionAction.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "math/GeometryComputations.hpp"
#include "utils/MasterSlave.hpp"

namespace precice {
namespace action {

logging::Logger BalanceVertexPositionAction::_log("action::BalanceVertexPositionAction");

BalanceVertexPositionAction:: BalanceVertexPositionAction
(
  Timing               timing,
  const mesh::PtrMesh& mesh,
  double               convergenceLimit,
  int                  maxIterations )
:
  Action(timing, mesh),
  _eps(convergenceLimit),
  _maxIterations(maxIterations)
{
  CHECK ( math::greater(_eps,0.0),
          "Convergence limit has to be larger than " << math::NUMERICAL_ZERO_DIFFERENCE
          << " for balance-vertex-position action!");
  CHECK ( maxIterations > 0,
          "Maximum iteration number has to be larger than 0 for "
          << "balance-vertex-position action!" );
}


void BalanceVertexPositionAction:: performAction
(
  double time,
  double dt,
  double computedPartFullDt,
  double fullDt )
{
  TRACE(dt, computedPartFullDt, fullDt);
  CHECK(not utils::MasterSlave::_masterMode && not utils::MasterSlave::_slaveMode,
        "BalanceVertexPositionAction is not yet supported for a usage with a Master")
  mesh::PtrMesh mesh = getMesh();
  assertion(mesh.get() != nullptr);
  using Eigen::VectorXd;
  using Eigen::Vector2d;
  using Eigen::Vector3d;
  int dimension = mesh->getDimensions();
  double errorMeasure = _eps * 10.0;
  int iterations = 0;
  while ( (errorMeasure > _eps) && (iterations < _maxIterations) ){
    errorMeasure = 0.0;
    std::vector<VectorXd> pullVectors ( mesh->vertices().size(), VectorXd::Zero(dimension) );
    VectorXd pullVectorScalings = VectorXd::Zero(mesh->vertices().size());
    // Gather pull vectors
    if ( dimension == 2 ) {
      for (mesh::Edge& edge : mesh->edges() ) {
        mesh::Vertex& v0 = edge.vertex(0);
        mesh::Vertex& v1 = edge.vertex(1);
        Vector2d ab = v0.getCoords();
        ab -= v1.getCoords();
        assertion ( v0.getID() < pullVectors.size(), v0.getID() );
        assertion ( v1.getID() < pullVectors.size(), v1.getID() );
        pullVectors[v0.getID()] -= ab;
        pullVectors[v1.getID()] += ab;
        pullVectorScalings[v0.getID()] += 1.0;
        pullVectorScalings[v1.getID()] += 1.0;
      }
    }
    else {
      assertion ( dimension == 3 );
      for (const mesh::Triangle& triangle : mesh->triangles()) {
        const mesh::Vertex& v0 = triangle.vertex(0);
        const mesh::Vertex& v1 = triangle.vertex(1);
        const mesh::Vertex& v2 = triangle.vertex(2);
        assertion ( v0.getID() < pullVectors.size(), v0.getID() );
        assertion ( v1.getID() < pullVectors.size(), v1.getID() );
        assertion ( v2.getID() < pullVectors.size(), v2.getID() );
        double area = math::GeometryComputations::triangleArea (
          v0.getCoords(),
          v1.getCoords(),
          v2.getCoords());
        Vector3d pullVector = triangle.getCenter();
        pullVector -= v0.getCoords();
        pullVector *= area;
        pullVectors[v0.getID()] += pullVector;
        pullVector = triangle.getCenter();
        pullVector -= v1.getCoords();
        pullVector *= area;
        pullVectors[v1.getID()] += pullVector;
        pullVector = triangle.getCenter();
        pullVector -= v2.getCoords();
        pullVector *= area;
        pullVectors[v2.getID()] += pullVector;
        pullVectorScalings[v0.getID()] += area;
        pullVectorScalings[v1.getID()] += area;
        pullVectorScalings[v2.getID()] += area;
      }
    }
    // Scale pull vectors
    for ( size_t i=0; i < mesh->vertices().size(); i++ ){
      pullVectors[i] /= pullVectorScalings[i];
    }
    // Project pull vector to get vertex displacement
    for (mesh::Vertex& vertex : mesh->vertices()) {
      if ( dimension == 2 ){
        Vector2d point = vertex.getCoords();
        assertion ( vertex.getID() < pullVectors.size(), vertex.getID() );
        point += pullVectors[vertex.getID()];
        Vector2d a = vertex.getCoords();
        // Get parameters for parametric representation of line orthogonal to vertex
        // normal: p(s) = a + s(b-a)
        Vector2d ab ( -vertex.getNormal()[1], vertex.getNormal()[0] );
        Vector2d b = a;
        b += ab;
        // Same for intersecting normal from point to line: q(t) = c + t(d - c)
        Vector2d c = point;
        Vector2d d = point;
        d += vertex.getNormal();
        // Compute denominator for solving 2x2 equation system
        double D = a(0)*(d(1)-c(1)) + b(0)*(c(1)-d(1)) + d(0)*ab(1) - c(0)*ab(1);
        assertion ( not math::equals(D, 0.0), a, b, c, d, ab );   // D == 0 would imply "normal // edge"
        // Compute parameter of intersection point on line ab
        double param = (a(0)*(d(1)-c(1)) + c(0)*(a(1)-d(1)) + d(0)*(c(1)-a(1))) / D;

        // Compute coordinates of projected point:
        VectorXd projectedPoint = ab;
        projectedPoint *= param;
        projectedPoint += a;

        // Compute error measure
        Vector2d coordDelta ( projectedPoint );
        coordDelta -= vertex.getCoords();
        errorMeasure += coordDelta.dot(coordDelta);
        vertex.setCoords ( projectedPoint );
      }
      else {
        assertion ( dimension == 3, dimension );
        Vector3d point = vertex.getCoords();
        assertion ( vertex.getID() < pullVectors.size(), vertex.getID() );
        point += pullVectors[vertex.getID()];
        Vector3d a = vertex.getCoords();
        // Parametric representation for plane orthogonal to vertex normal:
        // (x, y, z) * normal = d
        Vector3d normal = vertex.getNormal ();
        double d = normal.dot(a);

        // Parametric description of line from point, orthogonal to plane:
        // point + t * normal = x     (where t is parameter)
        // Determine t such that x lies on plane:
        double t = d - point.dot(normal) / normal.dot(normal);

        // Compute projected point with parameter t:
        VectorXd projectedPoint = vertex.getNormal();
        projectedPoint *= t;
        projectedPoint += point;

        // Compute error measure
        Vector3d coordDelta = projectedPoint;
        coordDelta -= vertex.getCoords();
        errorMeasure += coordDelta.dot(coordDelta);
        vertex.setCoords ( projectedPoint );
      }
    }
    mesh->computeState();
    //INFO ( "Error measure = " << errorMeasure );
    errorMeasure = std::sqrt(errorMeasure);
    iterations ++;
  }
  mesh->notifyListeners();
}

}} // namespace precice, action
