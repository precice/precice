#pragma once

#include <Eigen/Core>
#include <iosfwd>
#include <string>
#include <vector>
#include "io/ExportXML.hpp"
#include "logging/Logger.hpp"

namespace precice::mesh {
class Mesh;
class Edge;
class Triangle;
} // namespace precice::mesh

namespace precice::io {

/** Exporter for VTU and PVTU.
 *
 * Writes meshes to VTU piece files.
 * Parallel participants additionally write a PVTU file.
 * The naming scheme allows to import these files into Paraview as time series.
 */
class ExportVTU : public ExportXML {
public:
  ExportVTU(
      std::string_view  participantName,
      std::string_view  location,
      const mesh::Mesh &mesh,
      ExportKind        kind,
      int               frequency,
      int               rank,
      int               size);

private:
  mutable logging::Logger _log{"io::ExportVTU"};

  std::string getVTKFormat() const override;
  std::string getParallelExtension() const override;
  std::string getPieceExtension() const override;
  std::string getPieceAttributes(const mesh::Mesh &mesh) const override;

  void writeParallelCells(std::ostream &out) const override;

  void exportConnectivity(std::ostream &outFile, const mesh::Mesh &mesh) const override;
};

} // namespace precice::io
