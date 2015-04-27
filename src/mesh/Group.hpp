// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_MESH_GROUP_HPP_
#define PRECICE_MESH_GROUP_HPP_

#include "mesh/SharedPointer.hpp"
#include "utils/Helpers.hpp"
#include "utils/PointerVector.hpp"

namespace precice {
  namespace mesh {
    class Vertex;
    class Edge;
    class Triangle;
    class Quad;
  }
}

// ---------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace mesh {

/**
 * @brief Holds references to Triangle, Edge, and Vertex objects.
 *
 * A Group is not responsible for creating or destroying it's content, it is
 * just a grouping of mesh elements.
 */
class Group
{
public:

   typedef utils::ptr_vector<Vertex>   VertexContainer;
   typedef utils::ptr_vector<Edge>     EdgeContainer;
   typedef utils::ptr_vector<Triangle> TriangleContainer;
   typedef utils::ptr_vector<Quad>     QuadContainer;

   /**
    * @brief Adds a Vertex object to the Group.
    */
   void add ( Vertex& vertex );
   void add ( Vertex* vertex );

   /**
    * @brief Adds an Edge object to the Group.
    */
   void add ( Edge& edge );
   void add ( Edge* edge );

   /**
    * @brief Adds a Triangle object to the Group.
    */
   void add ( Triangle& triangle );
   void add ( Triangle* triangle );

   /**
    * @brief Adds a Quad object to the Group.
    */
   void add ( Quad& quad );
   void add ( Quad* quad );

   /**
    * @brief Adds the mesh elements of another Group object to the Group.
    */
   void add ( Group& group );

   /**
    * @brief Returns container with pointers to contained Vertex objects.
    */
   VertexContainer& vertices();

   /**
    * @brief Returns const container with pointers to contained Vertex objects.
    */
   const VertexContainer& vertices() const;

   /**
    * @brief Returns container with pointers to contained Edge objects.
    */
   EdgeContainer& edges();

   /**
    * @brief Returns const container with pointers to contained Edge objects.
    */
   const EdgeContainer& edges() const;

   /**
    * @brief Returns container with pointers to contained Triangle objects.
    */
   TriangleContainer& triangles();

   /**
    * @brief Returns const container with pointers to contained Triangle objects.
    */
   const TriangleContainer& triangles() const;

   /**
    * @brief Returns container with pointers to contained Quad objects.
    */
   QuadContainer& quads();

   /**
    * @brief Returns const container with pointers to contained Quad objects.
    */
   const QuadContainer& quads() const;

   /**
    * @brief Returns true, if no objects are contained in the group.
    */
   bool empty() const;

   /**
    * @brief Returns the number of all elements in the group.
    */
   size_t size() const;

   /**
    * @brief Removes all references to elements in the group.
    */
   void clear();

private:

   // @brief Container holding pointers to contained Vertex objects.
   VertexContainer _vertices;

   // @brief Container holding pointers to contained Edge objects.
   EdgeContainer _edges;

   // @brief Container holding pointers to contained Triangle objects.
   TriangleContainer _triangles;

   // @brief Container holding pointers to contained Quad objects.
   QuadContainer _quads;
};

}} // namespace precice, mesh

#endif /* PRECICE_MESH_GROUP_HPP_ */
