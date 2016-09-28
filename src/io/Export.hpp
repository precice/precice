#ifndef PRECICE_IO_EXPORT_HPP_
#define PRECICE_IO_EXPORT_HPP_

#include "Constants.hpp"
#include <string>

namespace precice {
  namespace mesh {
    class Mesh;
  }
}

// ---------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace io {

/**
 * @brief Abstract base class of all classes exporting container data structures.
 */
class Export
{
public:

  /**
   * @brief Constructor.
   *
   * @param location [IN] Location of export, filepath, e.g.
   * @param name [IN] Name of the export, filename, e.g.
   * @param exportable [IN] Vistable/s to be exported.
   */
  //Export();

  /**
   * @brief Destructor.
   */
  virtual ~Export() {}

  /**
   * @brief Returns the export type ID.
   */
  virtual int getType() const =0;

  /**
   * @brief Does export. Has to be implemented in subclass.
   *
   * @param name [IN] Filename (without path).
   * @param name [IN] Location (path without filename).
   * @param mesh [IN] Mesh to be exported.
   */
  virtual void doExport (
    const std::string& name,
    const std::string& location,
    mesh::Mesh&        mesh ) =0;
};

}} // namespace precice, io

#endif /* PRECICE_IO_EXPORT_HPP_ */
