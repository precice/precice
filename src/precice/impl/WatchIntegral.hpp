#pragma once

#include <Eigen/Core>
#include <limits>
#include <string>
#include <vector>
#include "mesh/SharedPointer.hpp"
#include "io/TXTTableWriter.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace mesh {
class Vertex;
}
} // namespace precice

namespace precice {
namespace impl {

class WatchIntegral {
public:
  /**
   * @brief Constructor.
   *
   * @param[in] meshToWatch Mesh to be watched, can be empty on construction.
   */
  WatchIntegral(
      mesh::PtrMesh      meshToWatch,
      const std::string &exportFilename);

  const mesh::PtrMesh &mesh() const;

  const std::string &filename() const;

  /// Writes one line with data of the integral over the mesh into the output file.
  void exportIntegralData(double time);

private:
  logging::Logger _log{"impl::WatchIntegral"};
  
  mesh::PtrMesh _mesh;

  io::TXTTableWriter _txtWriter;

  std::vector<mesh::PtrData> _dataToExport;

  void calculateVectorData(Eigen::VectorXd& value, mesh::PtrData data);

  void calculateScalarData(double& value, mesh::PtrData data);

  double calculateSurfaceArea();
};

} // namespace impl
} // namespace precice
