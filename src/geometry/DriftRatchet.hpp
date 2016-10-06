#ifndef PRECICE_GEOMETRY_DRIFTRATCHET_HPP_
#define PRECICE_GEOMETRY_DRIFTRATCHET_HPP_

#include "Geometry.hpp"

namespace precice {
   namespace mesh {
      class Mesh;
      class PropertyContainer;
      class Vertex;
      class Edge;
      class Triangle;
   }
}

// ---------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace geometry {

/**
 * A Drift Ratchet
 *
 * Represents the built-in drift ratchet geometry from Schindler.
 * This geometry runs along the x-axis and is axially symmetric.
 * The position of the geometry is fixed as soon as the rotation axis is fixed.
 *
 * !!! Parametrisation
 *
 * The original formula for the drift ratchet has been
 *
 * @f$ R(z) = R_min + (R_max - R_min) * (1 + g(z/L, k)) / 2 @f$
 *
 * with
 *
 * @f$ g(t,k) = sin(2*pi*t - k*g(t,k)) @f$
 *
 * Typical values for the input parameters are
 *
|| parameter || value
|  R_min     |  1.25
|  R_max     |  2.40
|  L         |  8.40
|  k         |  0.61
 *
 * I provide operations that derive all the parameters from the maxradius. See
 * the class' static methods. The parameter k is called shapeParameter.
 *
 * @author Tobias Weinzierl
 * @version $Revision: 1.12 $
 */
class DriftRatchet : public Geometry
{
public:

  /**
   * Construct Drift Ratchet
   *
   * @param axis Defines the offset (starting point) of the axis of the
   *             axially symmetric shape.
   */
  DriftRatchet(
    const Eigen::VectorXd&  offset,
    double                  discretizationWidth,
    double                  maxRadius,
    double                  minRadius,
    double                  shapeParameter,
    double                  length,
    double                  pores,
    int                     wallIndex,
    int                     inflowIndex,
    int                     outflowIndex );

  virtual ~DriftRatchet ();

  /**
   * This method returns the hard coded characteristic length L from the
   * formula @f$ \frac{r_{\textrm{max}}}{2.40} \cdot 8.4 @f$.
   */
  static double getCharacteristicLength ( double maxRadius );

  static double getDefaultMinRadius ( double maxRadius );

  static double getDefaultMaxRadius ();

  static double getDefaultShapeParameter ();

protected:

  virtual void specializedCreate ( mesh::Mesh& seed );

private:

  /**
   * The radius at a given depth z has to be computed recursively. This
   * constant determines the number of recursive steps.
   */
  static const int MAX_RECURSION_DEPTH;
  static const int INFLOW_GEO_ID;
  static const int OUTFLOW_GEO_ID;
  static const int WALL_GEO_ID;

  double _discretizationWidth;
  double _maxRadius;
  double _minRadius;
  double _shapeParameter;
  double _length;
  double _pores;

  /**
   * This method is used to compute g (?).
   *
   * It is used in getRadius().
   */
  double getG (
    double normalisedX,
    int    remainingIterations ) const;

  /**
   * Calculate the Radius
   *
   * The operation is given a scaled x, i.e. the parameter equals z/L.
   * @return Radius.
   */
  double getRadius ( double normalisedX ) const;

  /**
   * Calculate Number of Vertices per Cut
   *
   * Get the maximum circumference at divide it by h.
   */
  int getNumberOfVerticesPerCut ( double h ) const;

  void createLeftWall (
    mesh::Mesh&              mesh,
    mesh::PropertyContainer* propertyContainer,
    mesh::Vertex*            cutVertices[],
    mesh::Edge*              cutEdges[]  );

  void createBodyWall (
    mesh::Mesh&              mesh,
    mesh::PropertyContainer* propertyContainer,
    mesh::Vertex*            cutVertices[],
    mesh::Edge*              cutEdges[] );

  void createBodyWallSection (
    mesh::Mesh&              mesh,
    mesh::PropertyContainer* propertyContainer,
    mesh::Vertex*            cutVertices[],
    mesh::Edge*              cutEdges[],
    const Eigen::VectorXd&   center,
    double                   radius);

  void createRightWall ( mesh::Mesh              & mesh,
    mesh::PropertyContainer* propertyContainer,
    mesh::Vertex*            cutVertices[],
    mesh::Edge*              cutEdges[]  );

  int getCutsAlongXAxis ( double discretisationWidth ) const;

  double getExtremeCoordinateInAxisDirection ( int n ) const;
};

}} // namespace precice, geometry

#endif // PRECICE_GEOMETRY_DRIFTRATCHET_HPP_
