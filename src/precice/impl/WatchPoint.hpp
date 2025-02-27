#pragma once

#include <Eigen/Core>
#include <limits>
#include <string>
#include <vector>
#include "SharedPointer.hpp"
#include "io/TXTTableWriter.hpp"
#include "logging/Logger.hpp"
#include "mapping/Polation.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice::mesh {
class Vertex;
}

namespace precice::impl {

/// Observes and exports coordinates of a point on the geometry.
class WatchPoint {
public:
  /**
   * @brief Constructor.
   *
   * @param[in] meshToWatch Mesh to be watched, can be empty on construction.
   */
  WatchPoint(
      Eigen::VectorXd    pointCoords,
      mesh::PtrMesh      meshToWatch,
      const std::string &exportFilename);

  const mesh::PtrMesh &mesh() const;

  const std::string &filename() const;

  /** Initializes the watch point for exporting point data.
   *
   * This can be called repeatedly to reinitialize the WatchPoint.
   */
  void initialize();

  /// Writes one line with data of the watchpoint into the output file.
  void exportPointData(double time);

  bool isClosest() const
  {
    return _isClosest;
  }

private:
  logging::Logger _log{"impl::WatchPoint"};

  Eigen::VectorXd _point;

  mesh::PtrMesh _mesh;

  io::TXTTableWriter _txtWriter;

  double _shortestDistance = std::numeric_limits<double>::max();

  std::unique_ptr<mapping::Polation> _interpolation;

  std::vector<mesh::PtrData> _dataToExport;

  /// Holds the information if this processor is the closest
  bool _isClosest = true;

  void getValue(
      Eigen::VectorXd &value,
      mesh::PtrData &  data,
      double           time);

  void getValue(
      double &       value,
      mesh::PtrData &data,
      double         time);
};

} // namespace precice::impl
