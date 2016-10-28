#pragma once

#include "Geometry.hpp"
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"
#include <Eigen/Dense>

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
     const Eigen::VectorXd& offset,
     double                 discretizationWidth,
     const Eigen::VectorXd& sidelengths );

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
  Eigen::VectorXd _sidelengths;
};

}} // namespace precice, geometry
