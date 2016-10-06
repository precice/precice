#ifndef PRECICE_GEOMETRY_SOLVERGEOMETRY_HPP_
#define PRECICE_GEOMETRY_SOLVERGEOMETRY_HPP_

#include "Geometry.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace geometry {

/**
 * @brief Geometry for those meshes that are not communicated (meshes with "provide"
 * but no counterpart "from"). This class is a simple wrapper without any functionality.
 */
class SolverGeometry : public Geometry
{
public:

  SolverGeometry (
    const Eigen::VectorXd& offset);

  virtual ~SolverGeometry() {}

protected:

  /**
   * @brief Is called from Geometry. Nothing to be done here.
   */
  virtual void specializedCreate ( mesh::Mesh& seed );

private:

  // @brief Logging device.
  static logging::Logger _log;
};

}} // namespace precice, geometry

#endif /* PRECICE_GEOMETRY_SOLVERGEOMETRY_HPP_ */
