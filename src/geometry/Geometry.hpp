#ifndef PRECICE_GEOMETRY_GEOMETRY_HPP_
#define PRECICE_GEOMETRY_GEOMETRY_HPP_

#include "logging/Logger.hpp"
#include <Eigen/Dense>

namespace precice {
   namespace mesh {
      class Mesh;
   }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace geometry {

/**
 * @brief Abstract base class for all geometries.
 *
 * A geometry is there to create geometrical objects, based on the data
 * structure defined in component container. It uses the Mesh class to setup
 * a specific kind of geometry. In addition, every geometry is enanced by a
 * unique name and id.
 */
class Geometry
{
public:

  /**
   * @brief Constructor.
   *
   * @param[in] name  Unique name for the geometry.
   * @param[in] isVolumeEnclosed  If true, the volume of the geometry is inside of its (closed) surface.
   * @param[in] offset Constant offset of all vertex coordinates.
   */
  Geometry ( const Eigen::VectorXd& offset );

  /**
   * @brief Destructor.
   */
  virtual ~Geometry() {}

  /**
   * @brief Prepares the geometry creation.
   * This function is only overwritten by CommunicatedGeometry
   */
  virtual void prepare ( mesh::Mesh& seed ){}

  /**
   * @brief Creates the geometry into the given seed Mesh.
   */
  void create ( mesh::Mesh& seed );

  /**
   * @brief Sets an offset from zero for the geometry.
   */
  void setOffset ( const Eigen::VectorXd& offset )
  {
    _offset = offset;
  }

  /**
   * @brief Returns the offset of the geometry from zero.
   */
  const Eigen::VectorXd& getOffset () const
  {
    return _offset;
  }


protected:

  /**
   * @brief Is called from create() to actually create the geometry.
   */
  virtual void specializedCreate ( mesh::Mesh& seed ) =0;

  /**
   * @brief By default, mesh data is allocated after specializedCreate.
   *
   * If a subclass of Geometry already allocated the mesh data and
   * assigned values, this method has to be overwritten to not to allocate and
   * overwrite the data values.
   */
  virtual void allocateDataValues ( mesh::Mesh& mesh );

private:

  static logging::Logger _log;

  // @brief Offset of reference point of geometry from zero point
  Eigen::VectorXd _offset;

};

}} // namespace precice, geometry

#endif /* PRECICE_GEOMETRY_GEOMETRY_HPP_ */
