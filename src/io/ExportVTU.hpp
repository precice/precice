#pragma once

#include <Eigen/Core>
#include <iosfwd>
#include <string>
#include <vector>
#include "io/ExportXML.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace mesh {
class Mesh;
class Edge;
class Triangle;
} // namespace mesh
} // namespace precice

namespace precice {
namespace io {

/** Exporter for VTU and PVTU.
 *
 * Writes meshes to VTU piece files.
 * Parallel participants additionally write a PVTU file.
 * The naming scheme allows to import these files into Paraview as time series.
 */
class ExportVTU : public ExportXML {
private:
  mutable logging::Logger _log{"io::ExportVTU"};

  std::string getVTKFormat() const override;
  std::string getMasterExtension() const override;
  std::string getPieceExtension() const override;
  std::string getPieceAttributes(const mesh::Mesh &mesh) const override;

  void writeMasterCells(std::ostream &out) const override;

  void exportConnectivity(std::ostream &outFile, const mesh::Mesh &mesh) const override;
};

} // namespace io
} // namespace precice
