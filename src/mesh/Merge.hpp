#pragma once

#include "mesh/Group.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Quad.hpp"
#include <set>

// ---------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace mesh {

/**
 * @brief Merges all visitable objects into one Group while avoiding duplicates.
 *
 * Works for Triangle, Edge, and Vertex objects at the moment.
 *
 * The uniqueness of the merge is ensured by creating and checking a
 * property for every merged visitable. If the same visitable is encountered
 * a second time, it already has the property and is not merged again, hence.
 * The created properties are deleted on destruction of the MergeUniquelyVisitor
 * object.
 */
class Merge
{
public:

  /**
   * @brief Constructor.
   */
  Merge();

  /**
   * @brief Destructor. Removes temporary property from merged objects.
   */
  ~Merge();

  /**
   * @brief Merges the content of a CONTAINER_T with already merged content.
   */
  template<typename CONTAINER_T>
  Group& operator() ( CONTAINER_T& container );

  /**
   * @brief Returns the merged visitables
   */
  Group& content();

private:

  // @brief Merged visitables
  Group _merged;
};

// --------------------------------------------------------- HEADER DEFINITIONS

template<typename CONTAINER_T>
Group& Merge:: operator() ( CONTAINER_T& container )
{
  _merged.clear();
  std::set<Vertex*> vertices;
  for ( Vertex& vertex : container.vertices() ){
    if ( vertices.find(&vertex) == vertices.end() ){
      vertices.insert(&vertex);
      _merged.add(vertex);
    }
  }
  vertices.clear();
  std::set<Edge*> edges;
  for ( Edge& edge : container.edges() ){
    if ( edges.find(&edge) == edges.end() ){
      edges.insert(&edge);
      _merged.add(edge);
    }
  }
  edges.clear();
  std::set<Triangle*> triangles;
  for ( Triangle& triangle : container.triangles() ){
    if ( triangles.find(&triangle) == triangles.end() ){
      triangles.insert(&triangle);
      _merged.add(triangle);
    }
  }
  triangles.clear();
  std::set<Quad*> quads;
  for ( Quad& quad : container.quads() ){
    if ( quads.find(&quad) == quads.end() ){
      quads.insert(&quad);
      _merged.add(quad);
    }
  }
  quads.clear();
  return _merged;
}

}} // namespace precice, mesh

