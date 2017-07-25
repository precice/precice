#pragma once

#include "Export.hpp"
#include "logging/Logger.hpp"
#include <string>
#include <Eigen/Core>

namespace precice {
   namespace mesh {
      class Mesh;
   }
}

namespace precice {
namespace io {

/// Writes polygonal, triangle, or quadrangle meshes to vtk files.
class ExportVTK : public Export
{
public:

  /**
   * @brief Standard constructor
   *
   * @param[in] filename  Name of the vtk file (including file extension)
   * @param[in] container Container holding geometry to be visualized
   */
  explicit ExportVTK ( bool exportNormals );

  /// Returns the VTK type ID.
  virtual int getType() const;

  /// Perform writing to vtk file
  virtual void doExport (
    const std::string& name,
    const std::string& location,
    mesh::Mesh&        mesh );

  static void initializeWriting (
    std::ofstream&     filestream );

  static void writeHeader ( std::ostream& outFile );

  static void writeVertex (
    const Eigen::VectorXd&  position,
    std::ostream&           outFile );

  static void writeLine (
    int           vertexIndices[2],
    std::ostream& outFile );

  static void writeTriangle (
    int           vertexIndices[3],
    std::ostream& outFile );

  static void writeQuadrangle (
    int           vertexIndices[4],
    std::ostream& outFile );

private:

   // @brief Logging device.
   static logging::Logger _log;

   // @brief By default set true: plot vertex normals, false: no normals plotting
   bool _writeNormals;

   void openFile (
    std::ofstream&     outFile,
    const std::string& filename ) const;

   void exportGeometry (
     std::ofstream& outFile,
     mesh::Mesh&    mesh );

   void exportData (
     std::ofstream& outFile,
     mesh::Mesh&    mesh );
};

}} // namespace precice, io
