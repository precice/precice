#ifndef PRECICE_QUERY_FINDCLOSESTQUAD_HPP_
#define PRECICE_QUERY_FINDCLOSESTQUAD_HPP_

#include "utils/Globals.hpp"
#include <array>
#include <limits>

namespace precice {
  namespace mesh {
    class Mesh;
    class Quad;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace query {

/**
 * @brief Finds the closest Quad object contained in a Mesh object.
 */
class FindClosestQuad
{
public:

  /**
   * @brief Constructor.
   *
   * @param[in] searchPoint Origin from which the closest Quad object should be found.
   */
  FindClosestQuad ( const Eigen::VectorXd& searchPoint );

  /**
   * @brief Searches for the closest quad on the given mesh object.
   *
   * @return True, if a closest quad could be found.
   */
  template<typename CONTAINER_T>
  bool operator() ( CONTAINER_T & container );

  /**
   * @brief Returns the coordinates of the search point.
   */
  const Eigen::VectorXd& getSearchPoint() const;

  /**
   * @brief Returns true, if a closest quad has been found.
   */
  bool hasFound() const;

  /**
   * @brief Returns the distance to the found quad.
   *
   * Precondition: Find has been called and returned true.
   */
  double getEuclidianDistance();

  /**
   * @brief Returns the found quad.
   *
   * Precondition: Find has been called and returned true.
   */
  mesh::Quad& getClosestQuad();

  /**
   * @brief Returns the vector from the search point to the projection point.
   *
   * Precondition: Find has been called and returned true.
   */
  const Eigen::VectorXd& getVectorToProjectionPoint() const;

  /**
   * @brief Returns parametric description value (index 0, 1, 2) of proj. point.
   */
  double getProjectionPointParameter ( int index ) const;

private:

  static logging::Logger _log;

  /// Search point coordinates.
  Eigen::VectorXd _searchPoint;

  // @brief Shortest distance to the found Quad object.
  double _shortestDistance;

  // @brief Vector from search point to projection point.
  Eigen::VectorXd _vectorToProjectionPoint;

  // @brief Quad coordinates of the projection point.
  std::array<double,4> _parametersProjectionPoint; // Does this make sense?

  // @brief Pointer to found Quad object.
  mesh::Quad* _closestQuad;

  void find ( mesh::Quad& quad );
};

// --------------------------------------------------------- HEADER DEFINITIONS

template<typename CONTAINER_T>
bool FindClosestQuad:: operator() ( CONTAINER_T& container )
{
  for ( mesh::Quad& quad : container.quads() ) {
    find ( quad );
  }
  return _closestQuad != NULL;
}

}} // namespace precice, query

#endif /* PRECICE_QUERY_FINDCLOSESTQUAD_HPP_ */

