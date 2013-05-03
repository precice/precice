// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ImportGeometry.hpp"
#include "io/Import.hpp"
#include "io/ImportVRML.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/Globals.hpp"

namespace precice {
namespace geometry {

ImportGeometry:: ImportGeometry
(
  const utils::DynVector& offset,
  const std::string&      fileName,
  FileType                fileType,
  bool                    importCheckpoint )
:
   Geometry ( offset ),
   _fileName ( fileName ),
   _fileType ( fileType ),
   _importCheckpoint ( importCheckpoint )
{}

void ImportGeometry:: specializedCreate
(
  mesh::Mesh& seed )
{
  if ( _fileType == VRML_1_FILE ){
    io::ImportVRML import ( "" );
    if ( _importCheckpoint ){
      import.doImportCheckpoint ( _fileName, seed );
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
