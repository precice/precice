#pragma once

#include <Eigen/Core>
#include <iosfwd>
#include <string>
#include <vector>
#include "io/Export.hpp"
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {
namespace mesh {
class Mesh;
class Edge;
class Triangle;
class Tetrahedron;
} // namespace mesh
} // namespace precice

namespace precice {
namespace io {

/// Common class to generate the VTK XML-based formats.
class ExportXML : public Export {
public:
  ExportXML(
      std::string_view  participantName,
      std::string_view  location,
      const mesh::Mesh &mesh,
      ExportKind        kind,
      int               frequency,
      int               rank,
      int               size);

  void doExport(int index, double time) final override;

  void exportSeries() const final override;

  static void writeVertex(
      const Eigen::VectorXd &position,
      std::ostream &         outFile);

  static void writeLine(
      const mesh::Edge &edge,
      std::ostream &    outFile);

  static void writeTriangle(
      const mesh::Triangle &triangle,
      std::ostream &        outFile);

  static void writeTetrahedron(
      const mesh::Tetrahedron &tetra,
      std::ostream &           outFile);

private:
  mutable logging::Logger _log{"io::ExportXML"};

  /// List of names of all scalar data on mesh
  std::vector<std::string> _scalarDataNames;

  /// List of names of all vector data on mesh
  std::vector<std::string> _vectorDataNames;

  virtual std::string getVTKFormat() const                             = 0;
  virtual std::string getParallelExtension() const                     = 0;
  virtual std::string getPieceExtension() const                        = 0;
  virtual std::string getPieceAttributes(const mesh::Mesh &mesh) const = 0;

  /**
   * @brief Stores scalar and vector data names in string vectors
   * Needed for writing primary file and sub files
   */
  void processDataNamesAndDimensions(const mesh::Mesh &mesh);

  /**
   * @brief Writes the primary file (called only by the primary rank)
   */
  void writeParallelFile(int index, double time);

  virtual void writeParallelCells(std::ostream &out) const = 0;

  void writeParallelData(std::ostream &out) const;

  /**
   * @brief Writes the sub file for each rank
   */
  void writeSubFile(int index, double time);

  void exportPoints(
      std::ostream &    outFile,
      const mesh::Mesh &mesh) const;

  virtual void exportConnectivity(
      std::ostream &    outFile,
      const mesh::Mesh &mesh) const = 0;

  void exportData(
      std::ostream &    outFile,
      const mesh::Mesh &mesh) const;

  void exportGradient(const mesh::PtrData data, const int dataDim, std::ostream &outFile) const;

  std::string parallelPieceFilenameFor(int index, int rank) const;
  std::string serialPieceFilename(int index) const;
};

} // namespace io
} // namespace precice
