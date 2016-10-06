#ifndef PRECICE_GEOMETRY_BUBBLE_HPP_
#define PRECICE_GEOMETRY_BUBBLE_HPP_

#include "Geometry.hpp"
#include "logging/Logger.hpp"
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
 * @brief Creates a triangulated mesh with an ellipse like shape (2D only).
 * This corresponds to the limit shape of an oscillating bubble.
 *
 * The elements to create the mesh are Mesh, Triangle and Vertex objects.
 *
 * The shape formula is taken from the bubble oscillation tests in
 *   A novel representation of the surface tension force for two-phase flow
 *     with reduced spurious currents
 *   by: E. Aulisa, S. Manservisi, R. Scardovelli
 *   in: Computer Methods in Applied Mechanics and Engineering, Vol. 195, No. 44-47. (2006)
 *
 * @author Florian Ferstl, Bernhard Gatzhammer
 *
 */
class Bubble : public Geometry
{
public:

   /**
    * @brief Standard constructor
    *
    * @param[in] radius      Radius of the (undeformed) bubble
    * @param[in] offset      Offset of the bubble's center from the origin
    * @param[in] deformation Deformation of bubble. Between -4/3 (= "8" like
    *                        shape) and 4/3 (= "oo" like shape). Value 0 results
    *                        in a perfect sphere with specified radius.
    */
  Bubble ( const Eigen::VectorXd&   offset,
            double                  discretizationWidth,
            double                  radius,
            double                  deformation );

   /**
    * @brief Destructor.
    */
   virtual ~Bubble() {}

protected:

   /**
    * @brief Creates the triangulated mesh of the bubblee into a container
    */
  virtual void specializedCreate ( mesh::Mesh& seed );

private:

   // @brief Logging device.
   static logging::Logger _log;

   // @brief Length of mesh elements the bubble is made of
   double _discretizationWidth;

   // @brief Radius of the (undeformed) bubble
   double _radius;

   // @brief Deformation of bubble
   double _deformation;

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

#endif /*PRECICE_GEOMETRY_BUBBLE_HPP_*/
