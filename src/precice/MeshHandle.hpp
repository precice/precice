#pragma once

#include <vector>
#include <cstddef>

namespace precice {
  namespace mesh {
    class Mesh;
    class Group;
  }
  namespace impl {
    struct VertexIteratorImplementation;
    struct EdgeIteratorImplementation;
    struct TriangleIteratorImplementation;
    // no quad iterator yet
  }
}


// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {

/**
 * @brief Iterator over vertex coordinates and IDs.
 */
class VertexIterator
{
public:

  VertexIterator (
    const mesh::Group& content,
    bool               begin );

  VertexIterator ( const VertexIterator& toCopy );

  VertexIterator& operator= ( const VertexIterator& toAssign );

  ~VertexIterator();

  /// Postfix operator
  VertexIterator& operator++(int unused);

  /// Prefix operator
  VertexIterator& operator++();

  VertexIterator& operator*();

  int vertexID();

  const double* vertexCoords();

  bool operator!= ( const VertexIterator& vertexIterator );

private:

  impl::VertexIteratorImplementation* _impl;
};

/**
 * @brief Offers methods begin() and end() to iterate over all vertices.
 */
class VertexHandle
{
public:

  typedef VertexIterator const_iterator;

  /**
   * @brief Constructor, reference to mesh object holding vertices required.
   */
  VertexHandle ( const mesh::Group& content );

  /**
   * @brief Returns iterator to begin of the geometry's vertices.
   */
  VertexIterator begin() const;

  /**
   * @brief Returns iterator to end of the geometry's vertices.
   */
  VertexIterator end() const;

  std::size_t size() const;

private:

  // @brief Group instance holding vertices.
  const mesh::Group& _content;
};

class EdgeIterator
{
public:

  EdgeIterator (
    const mesh::Group& mesh,
    bool              begin );

  ~EdgeIterator();

  EdgeIterator& operator++(int);

  const double* vertexCoords ( int vertexIndex );

  int vertexID ( int vertexIndex );

  bool operator!= ( const EdgeIterator& edgeIterator );

private:

  impl::EdgeIteratorImplementation* _impl;
};

/**
 * @brief Offers methods begin() and end() to iterate over all edges.
 */
class EdgeHandle
{
public:

   // @brief Necessary to be used by boost::foreach.
   typedef EdgeIterator const_iterator;

   /**
    * @brief Constructor, reference to mesh object holding edges required.
    */
   EdgeHandle ( const mesh::Group& mesh );

   /**
    * @brief Returns iterator to begin of the geometry's edges.
    */
   EdgeIterator begin() const;

   /**
    * @brief Returns iterator to end of the geometry's edges.
    */
   EdgeIterator end() const;

  std::size_t size() const;

private:

   // @brief Mesh instance holding edges.
   const mesh::Group& _content;
};

class TriangleIterator
{
public:

  TriangleIterator (
    const mesh::Group& content,
    bool               begin );

  ~TriangleIterator();

  TriangleIterator& operator++(int);

  const double* vertexCoords ( int vertexIndex );

  int vertexID ( int vertexIndex );

  bool operator!= ( const TriangleIterator& triangleIterator );

private:

  impl::TriangleIteratorImplementation* _impl;
};

/**
 * @brief Offers methods begin() and end() to iterate over all triangles.
 */
class TriangleHandle
{
public:

   typedef TriangleIterator const_iterator;

   /**
    * @brief Constructor, reference to mesh object holding triangles required.
    */
   TriangleHandle ( const mesh::Group& content );

   /**
    * @brief Returns iterator to begin of the geometry's triangles.
    */
   TriangleIterator begin() const;

   /**
    * @brief Returns iterator to end of the geometry's triangles.
    */
   TriangleIterator end() const;

  std::size_t size() const;

private:

   // @brief Mesh instance holding triangles.
   const mesh::Group& _content;
};

/**
 * @brief Allows to query vertices, edges, and triangles of a geometry.
 *
 * A geometry handle can be retrieved from the coupling interface by the method
 * precice::SolverInterfaceImpl::getMeshHandle().
 *
 * Access to vertices is done via the method precice::MeshHandle::vertices(),
 * which returns a MeshHandle::VertexHandle object, which offers the
 * methods begin() and end() to iterate over all vertices. Access to edges and
 * triangles is analogous.
 *
 * The iterators return const references to the respective objects on dereferenciation
 * by using operator*(). The operator -> is not supported, and does only return
 * pointers to the objects.
 */
class MeshHandle
{
public:

   /**
    * @brief Standard constructor, not meant to be used by a solver.
    *
    * @param mesh [IN] The mesh representing the geometry.
    */
   MeshHandle ( const mesh::Group& content );

   /**
    * @brief Returns handle for Vertex objects.
    */
   const VertexHandle& vertices () const;

   /**
    * @brief Returns handle for Edge objects.
    */
   const EdgeHandle& edges () const;

   /**
    * @brief Returns handle for Triangle objects.
    */
   const TriangleHandle& triangles () const;

private:

   // @brief Handle for vertices.
   VertexHandle _vertexHandle;

   // @brief Handle for edges.
   EdgeHandle _edgeHandle;

   // @brief Handle for triangles.
   TriangleHandle _triangleHandle;
};

} // namespace precice
