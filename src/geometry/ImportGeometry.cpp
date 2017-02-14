#include "ImportGeometry.hpp"
#include "io/Import.hpp"
#include "io/ImportVRML.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {
namespace geometry {

ImportGeometry:: ImportGeometry
(
  const Eigen::VectorXd&  offset,
  const std::string&      fileName,
  FileType                fileType,
  bool                    importCheckpoint,
  bool                    createMesh)
:
   Geometry ( offset ),
   _fileName ( fileName ),
   _fileType ( fileType ),
   _importCheckpoint ( importCheckpoint ),
   _createMesh( createMesh )
{}

void ImportGeometry:: specializedCreate
(
  mesh::Mesh& seed )
{
  if ( _fileType == VRML_1_FILE ){
    io::ImportVRML import ( "" );
    if ( _importCheckpoint ){
      import.doImportCheckpoint ( _fileName, seed, _createMesh );
    }
    else {
      import.doImport ( _fileName, seed );
    }
  }
}

void ImportGeometry:: allocateDataValues
(
  mesh::Mesh& mesh )
{
  // Data values are already allocated during import when a checkpoint is read
  // and should not be overwritten by reallocating them.
  if ( not _importCheckpoint ){
    mesh.allocateDataValues();
  }
}

}} // namespace precice, geometry
