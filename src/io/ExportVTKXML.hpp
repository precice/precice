#pragma once

#include "Export.hpp"
#include "logging/Logger.hpp"
#include <vector>
#include <string>
#include <Eigen/Dense>

namespace precice {
   namespace mesh {
      class Mesh;
      class Quad;
      class Edge;
      class Triangle;
   }
}

namespace precice {
namespace io {

/**
 * @brief Writes meshes to xml-vtk files. Only for parallel usage. Serial usage (coupling mode) should still use ExportVTK
 */
class ExportVTKXML : public Export
{
public:

  /**
   * @brief Standard constructor
   *
   * @param exportNormals  [IN] boolean: write normals to file?
   */
  ExportVTKXML ( bool writeNormals);

  /// Returns the VTK type ID.
  virtual int getType() const;

  /// Perform writing to vtk file
  virtual void doExport (
    const std::string& name,
    const std::string& location,
    mesh::Mesh&        mesh );

  static void writeVertex (
    const Eigen::VectorXd& position,
    std::ofstream&         outFile );

  static void writeLine (
    mesh::Edge&    edge,
    std::ofstream& outFile );

  static void writeTriangle (
    mesh::Triangle& triangle,
    std::ofstream&  outFile );

  static void writeQuadrangle (
    mesh::Quad&    quad,
    std::ofstream& outFile );

private:

   static logging::Logger _log;

   /// By default set true: plot vertex normals, false: no normals plotting
   bool _writeNormals;

   /// dimensions of mesh
   int _meshDimensions;

   /// List of names of all scalar data on mesh
   std::vector<std::string> _scalarDataNames;

   /// List of names of all vector data on mesh
   std::vector<std::string> _vectorDataNames;

   /**
    * @brief Stores scalar and vector data names in string vectors
    * Needed for writing master file and sub files
    */
   void processDataNamesAndDimensions
   (
     mesh::Mesh& mesh);

   /**
    * @brief Writes the master file (called only by the master rank)
    */
   void writeMasterFile
   (
     const std::string& name,
     const std::string& location,
     mesh::Mesh&        mesh);

   /**
    * @brief Writes the sub file for each rank
    */
   void writeSubFile
   (
     const std::string& name,
     const std::string& location,
     mesh::Mesh&        mesh);

   void exportGeometry (
     std::ofstream& outFile,
     mesh::Mesh&    mesh );

   void exportData (
     std::ofstream& outFile,
     mesh::Mesh&    mesh );
};

}} // namespace precice, io
