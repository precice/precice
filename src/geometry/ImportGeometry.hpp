#ifndef PRECICE_GEOMETRY_IMPORTGEOMETRY_HPP_
#define PRECICE_GEOMETRY_IMPORTGEOMETRY_HPP_

#include "Geometry.hpp"

namespace precice {
namespace geometry {

/**
 * @brief Imports a geometry from a file.
 */
class ImportGeometry
:
   public Geometry
{
public:

  /**
   * @brief Possible file types to be imported.
   */
  enum FileType {
   VRML_1_FILE
  };

  /**
   * @brief Constructor.
   *
   * @param name [IN] Unique name for the geometry.
   * @param isVolumeEnclosed [IN] Default normal direction, if true.
   * @param offset [IN] Homogeneous offset of mesh nodes.
   * @param fileName [IN] Path and name of the file to be imported.
   * @param fileType [IN] Type of the file to be imported.
   */
  ImportGeometry (
   const utils::DynVector& offset,
   const std::string&      fileName,
   FileType                fileType,
   bool                    importCheckpoint,
   bool                    createMesh);

  /**
   * @brief Destructor.
   */
  virtual ~ImportGeometry() {}

protected:

  virtual void specializedCreate ( mesh::Mesh& seed );

  virtual void allocateDataValues ( mesh::Mesh& mesh );

private:

   // @brief Locates the file to import with directory and name.
   std::string _fileName;

   // @brief Determines the type of the file to be imported.
   FileType _fileType;

   bool _importCheckpoint;
   bool _createMesh;
};

}} // namespace precice, geometry

#endif /* PRECICE_GEOMETRY_IMPORTGEOMETRY_HPP_ */
