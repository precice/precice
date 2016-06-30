#ifndef PRECICE_GEOMETRY_CUBOID_HPP_
#define PRECICE_GEOMETRY_CUBOID_HPP_

#include "Geometry.hpp"
#include "logging/Logger.hpp"
#include "utils/Dimensions.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {
namespace geometry {


class Cuboid : public Geometry
{
public:

   /**
    * @brief Constructor
    *
    * @param name        Name of the cuboid geometry
    * @param offset      Offset of corner with coordinates x(i) = 0, for
    *                    all i=1,..,Def::DIM from origin
    * @param sidelengths Lengths of the sides of the cuboid
    */
   Cuboid (
     const utils::DynVector& offset,
     double                  discretizationWidth,
     const utils::DynVector& sidelengths );

   /**
    * @brief Destructor
    */
   virtual ~Cuboid () {};

protected:

   virtual void specializedCreate ( mesh::Mesh& seed );

private:

   // @brief Logging device.
   static logging::Logger _log;

   // @brief Determines minimal length of elements used to build the Cuboid
   double _discretizationWidth;

   // @brief Sidelengths of the cuboid
   utils::DynVector _sidelengths;
};

}} // namespace precice, geometry

#endif /* PRECICE_GEOMETRY_CUBOID_HPP_ */
