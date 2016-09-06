#ifndef PRECICE_IO_EXPORTVTK_HPP_
#define PRECICE_IO_EXPORTVTK_HPP_

#include "Export.hpp"
#include "logging/Logger.hpp"
#include "tarch/la/Vector.h"
#include "utils/Dimensions.hpp"
#include <string>

namespace precice {
   namespace mesh {
      class Mesh;
   }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace io {

/**
 * @brief Writes polygonal, triangle, or quadrangle meshes to vtk files.
 */
class ExportVTK : public Export
{
public:

  /**
   * @brief Standard constructor
   *
   * @param filename  [IN] Name of the vtk file (including file extension)
   * @param container [IN] Container holding geometry to be visualized
   */
  ExportVTK ( bool exportNormals );

  /**
   * @brief Returns the VTK type ID.
   */
  virtual int getType() const;

  /**
   * @brief Perform writing to vtk file
   */
  virtual void doExport (
    const std::string& name,
    const std::string& location,
    mesh::Mesh&        mesh );

  static void initializeWriting (
    const std::string& filename,
    std::ofstream&     filestream );

  static void writeHeader ( std::ostream& outFile );

  static void writeVertex (
    const utils::DynVector& position,
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

#endif /* PRECICE_IO_EXPORTVTK_HPP_ */
