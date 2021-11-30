#pragma once

#include <Eigen/Core>
#include <iosfwd>
#include <string>
#include <vector>
#include "io/Export.hpp"
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

/// Common class to generate the VTK XML-based formats.
class ExportXML : public Export {
public:
  void doExport(
      const std::string &name,
      const std::string &location,
      const mesh::Mesh & mesh) override;

  static void writeVertex(
      const Eigen::VectorXd &position,
      std::ostream &         outFile);

  static void writeLine(
      const mesh::Edge &edge,
      std::ostream &    outFile);

  static void writeTriangle(
      const mesh::Triangle &triangle,
      std::ostream &        outFile);

private:
  mutable logging::Logger _log{"io::ExportXML"};

  /// List of names of all scalar data on mesh
  std::vector<std::string> _scalarDataNames;

  /// List of names of all vector data on mesh
  std::vector<std::string> _vectorDataNames;

  virtual std::string getVTKFormat() const                             = 0;
  virtual std::string getMasterExtension() const                       = 0;
  virtual std::string getPieceExtension() const                        = 0;
  virtual std::string getPieceAttributes(const mesh::Mesh &mesh) const = 0;

  /**
    * @brief Stores scalar and vector data names in string vectors
    * Needed for writing master file and sub files
    */
  void processDataNamesAndDimensions(const mesh::Mesh &mesh);

  /**
    * @brief Writes the master file (called only by the master rank)
    */
  void writeMasterFile(
      const std::string &name,
      const std::string &location,
      const mesh::Mesh & mesh) const;

  virtual void writeMasterCells(std::ostream &out) const = 0;

  void writeMasterData(std::ostream &out) const;

  /**
    * @brief Writes the sub file for each rank
    */
  void writeSubFile(
      const std::string &name,
      const std::string &location,
      const mesh::Mesh & mesh) const;

  void exportPoints(
      std::ostream &    outFile,
      const mesh::Mesh &mesh) const;

  virtual void exportConnectivity(
      std::ostream &    outFile,
      const mesh::Mesh &mesh) const = 0;

  void exportData(
      std::ostream &    outFile,
      const mesh::Mesh &mesh) const;
};

} // namespace io
} // namespace precice
