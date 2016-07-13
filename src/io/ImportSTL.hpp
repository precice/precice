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
