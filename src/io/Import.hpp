#ifndef PRECICE_IO_IMPORT_HPP_
#define PRECICE_IO_IMPORT_HPP_

#include "utils/Dimensions.hpp"
#include <string>
#include <map>

namespace precice {
  namespace mesh {
    class Mesh;
  }
}

// ------------------------------------------------------------ CLASS DEFINTION

namespace precice {
namespace io {

/**
 * @brief Abstract base class for import classes.
 *
 * Defines a unique interface to import geometries.
 */
class Import
{
public:

  Import ( const std::string& location );

  virtual ~Import() {}

  /**
   * @brief Does the import.
   *
   * @param mesh [IN/OUT] The importet elements are added to the mesh.
   */
  virtual void doImport (
    const std::string& name,
    mesh::Mesh&        mesh ) =0;

protected:

  const std::string& getLocation () const
  {
    return _location;
  }

private:

  std::string _location;
};

}} // namespace precice, io

#endif /* PRECICE_IO_IMPORT_HPP_ */
