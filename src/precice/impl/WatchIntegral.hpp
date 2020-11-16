#pragma once

#include <Eigen/Core>
#include <string>
#include <vector>
#include "io/TXTTableWriter.hpp"
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {
namespace mesh {
class Vertex;
}
} // namespace precice

namespace precice {
namespace impl {

/**
 * @brief Track and output transient integral data on a mesh
 * 
 * Calculation of integral depends on connectivity information of mesh and
 * given scaling option. If scale with area option is true, vertex data is 
 * weighted by area and summed up, otherwise, vertex data is directly summed
 * up.
 */
class WatchIntegral {
public:
  /**
   * @brief Constructor.
   *
   * @param[in] meshToWatch Mesh to be watched.
   * @param[in] exportFilename output file name
   * @param[in] isScalingOn whether the data will be scaled with area or not
   */
  WatchIntegral(
      mesh::PtrMesh      meshToWatch,
      const std::string &exportFilename,
      bool               isScalingOn);

  /// Writes one line with data of the integral over the mesh into the output file.
  void exportIntegralData(double time);

  /// Adds surface area information based on mesh connectivity
  void initialize();

private:
  logging::Logger _log{"impl::WatchIntegral"};

  mesh::PtrMesh _mesh;

  io::TXTTableWriter _txtWriter;

  std::vector<mesh::PtrData> _dataToExport;

  bool _isScalingOn;

  Eigen::VectorXd calculateVectorData(mesh::PtrData data);

  double calculateScalarData(mesh::PtrData data);

  double calculateSurfaceArea();
};

} // namespace impl
} // namespace precice
