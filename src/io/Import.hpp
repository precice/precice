#ifndef PRECICE_IO_IMPORT_HPP_
#define PRECICE_IO_IMPORT_HPP_

#include <map>
#include <string>

namespace precice {
namespace mesh {
class Mesh;
}
} // namespace precice

// ------------------------------------------------------------ CLASS DEFINTION

namespace precice {
namespace io {

/**
 * @brief Abstract base class for import classes.
 *
 * Defines a unique interface to import geometries.
 */
class Import {
public:
  Import(const std::string &location);

  virtual ~Import() {}

  /**
   * @brief Does the import.
   *
   * @param mesh [IN/OUT] The importet elements are added to the mesh.
   */
  virtual void doImport(
      const std::string &name,
      mesh::Mesh &       mesh) = 0;

protected:
  const std::string &getLocation() const
  {
    return _location;
  }

private:
  std::string _location;
};

} // namespace io
} // namespace precice

#endif /* PRECICE_IO_IMPORT_HPP_ */
