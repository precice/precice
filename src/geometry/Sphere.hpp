#ifndef PRECICE_GEOMETRY_SPHERE_HPP_
#define PRECICE_GEOMETRY_SPHERE_HPP_

#include "Geometry.hpp"
#include "utils/Helpers.hpp"
#include "mesh/SharedPointer.hpp"
#include "tarch/logging/Log.h"
#include <map>

namespace precice {
  namespace mesh {
    class Mesh;
    class Edge;
    class Vertex;
  }
}

// ---------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace geometry {


/**
 * @brief Creates a sphere (in 3D) or circle (in 2D) triangulated mesh
 *
 * The elements to create the mesh are Mesh, Triangle and Vertex objects.
 *
 * @author Bernhard Gatzhammer
 *
 */
class Sphere : public Geometry
{
public:

   /**
    * @brief Standard constructor
    *
    * @param radius  [IN] Radius of the sphere
    * @param offset  [IN] Offset of the sphere's center from the origin
    */
   Sphere (
     const utils::DynVector& offset,
     double                  discretizationWidth,
     double                  radius );

   /**
    * @brief Destructor.
    */
   virtual ~Sphere() {}

protected:

   /**
    * @brief Creates the triangulated mesh of the sphere into a container
    */
//   virtual mesh::PtrMesh specializedCreate ();

   virtual void specializedCreate ( mesh::Mesh& seed );

private:

   // @brief Logging device.
   static tarch::logging::Log _log;

   // @brief Length of mesh elements the sphere is made of
   double _discretizationWidth;

   // @brief Radius of the sphere
   double _radius;

   mesh::Vertex* getVertex (
     mesh::Vertex&                               v0,
     mesh::Vertex&                               v1,
     std::map<std::pair<int,int>,mesh::Vertex*>& dividedEdges,
     mesh::Mesh&                                 seed );

   mesh::Edge* getEdge (
     mesh::Vertex&                             v0,
     mesh::Vertex&                             v1,
     std::map<std::pair<int,int>,mesh::Edge*>& edges,
     mesh::Mesh&                               seed );

   int getNumberLongitudinalElements ( double discretizationWidth ) const;

   int getNumberLatitudinalElements ( double discretizationWidth ) const;
};

}} // namespace precice, geometry

#endif /*PRECICE_GEOMETRY_SPHERE_HPP_*/
