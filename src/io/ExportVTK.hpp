#pragma once

#include <Eigen/Core>
#include <iosfwd>
#include <string>
#include "Export.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace mesh {
class Mesh;
}
} // namespace precice

namespace precice {
namespace io {

/// Writes polygonal, or triangle meshes to vtk files.
class ExportVTK : public Export {
public:
  explicit ExportVTK(bool exportNormals);

  /// Returns the VTK type ID.
  virtual int getType() const;

  /// Perform writing to VTK file
  virtual void doExport(
      const std::string &name,
      const std::string &location,
      mesh::Mesh &       mesh);

  static void initializeWriting(
      std::ofstream &filestream);

  static void writeHeader(std::ostream &outFile);

  static void writeVertex(
      const Eigen::VectorXd &position,
      std::ostream &         outFile);

  static void writeLine(
      int           vertexIndices[2],
      std::ostream &outFile);

  static void writeTriangle(
      int           vertexIndices[3],
      std::ostream &outFile);

private:
  logging::Logger _log{"io::ExportVTK"};

  /// By default set true: plot vertex normals, false: no normals plotting
  bool _writeNormals;

  void openFile(
      std::ofstream &    outFile,
      const std::string &filename) const;

  void exportMesh(std::ofstream &outFile, mesh::Mesh const &mesh);

  void exportData(std::ofstream &outFile, mesh::Mesh const &mesh);
};

} // namespace io
} // namespace precice
