#ifndef PRECICE_QUERY_FINDCLOSESTVERTEX_HPP_
#define PRECICE_QUERY_FINDCLOSESTVERTEX_HPP_

#include "mesh/Vertex.hpp"
#include <limits>

// ---------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace query {

/**
 * @brief Finds the closest Vertex object hold by Mesh objects.
 */
class FindClosestVertex
{
public:

  /**
   * @brief Constructor.
   *
   * @param[in] searchPoint Coordinates of origin of search for closest point.
   */
  FindClosestVertex ( const Eigen::VectorXd& searchPoint );

  /**
   * @brief Searches among all Vertex objects hold by the given Mesh object.
   *
   * When called for different meshes, the closest Vertex object in among all
   * the meshes is found.
   */
  template<typename CONTAINER_T>
  bool operator() ( CONTAINER_T& container );

  /**
   * @brief Returns the coordinates of the search point.
   */
  const Eigen::VectorXd& getSearchPoint() const;

  bool hasFound() const;

  /**
   * @brief Returns the euclidian distance to the closest vertex found.
   *
   * Precondition: find() has been called and returned true.
   */
  double getEuclidianDistance();

  /**
   * @brief Returns the closest vertex found.
   *
   * Precondition: find() has been called and returned true.
   */
  mesh::Vertex& getClosestVertex();

private:

  // @brief Origin of search for closest vertex.
  Eigen::VectorXd _searchPoint;

  // @brief Distance to closest Vertex object found.
  double _shortestDistance;

  // @brief Pointer to closest Vertex object found.
  mesh::Vertex* _closestVertex;
};

// --------------------------------------------------------- HEADER DEFINITIONS

template<typename CONTAINER_T>
bool FindClosestVertex:: operator()
(
  CONTAINER_T& container )
{
  Eigen::VectorXd vectorDistance = Eigen::VectorXd::Zero(_searchPoint.size());
  for ( mesh::Vertex& vertex : container.vertices() ) {
    assertion ( vertex.getDimensions() == _searchPoint.size(),
                vertex.getDimensions(), _searchPoint.size() );
    vectorDistance = vertex.getCoords();
    vectorDistance -= _searchPoint;
    double distance = vectorDistance.norm();
    if ( distance < _shortestDistance) {
      _shortestDistance = distance;
      _closestVertex = &vertex;
    }
  }
  return _closestVertex != NULL;
}

}} // namespace precice, query

#endif /* PRECICE_QUERY_FINDCLOSESTVERTEX_HPP_ */

