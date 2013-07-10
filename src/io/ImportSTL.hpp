// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IO_IMPORTSTL_HPP_
#define PRECICE_IO_IMPORTSTL_HPP_

#include "Import.hpp"
#include "tarch/logging/Log.h"
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
 * @brief Imports a geometry from a vrml file.
 */
class ImportSTL : public Import
{
public:

   /**
    * @brief Constructor.
    */
  ImportSTL(const std::string& location);


   /**
    * @brief Destructor, empty.
    */
   virtual ~ImportSTL() {}

   /**
    * @brief Imports the geometry from an STL file into a Mesh object.
    *
    * @param mesh [IN/OUT] The imported elements are added to this mesh.
    */
   virtual void doImport (
      const std::string& name,
      mesh::Mesh&        mesh );

private:

   // @brief Logging device.
   static tarch::logging::Log _log;
};

}} // namespace precice, io

#endif /* PRECICE_IO_IMPORTSTL_HPP_ */
