#ifndef PRECICE_IO_EXPORT_HPP_
#define PRECICE_IO_EXPORT_HPP_

#include "Constants.hpp"
#include <string>

namespace precice {
  namespace mesh {
    class Mesh;
  }
}

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
   * @param[in] location  Location of export, filepath, e.g.
   * @param[in] name Name of the export, filename, e.g.
   * @param[in] exportable Vistable/s to be exported.
   */
  //Export();

  virtual ~Export() {}

  /// Returns the export type ID.
  virtual int getType() const =0;

  /**
   * @brief Does export. Has to be implemented in subclass.
   *
   * @param[in] name Filename (without path).
   * @param[in] name Location (path without filename).
   * @param[in] mesh Mesh to be exported.
   */
  virtual void doExport (
    const std::string& name,
    const std::string& location,
    mesh::Mesh&        mesh ) =0;
};

}} // namespace precice, io

#endif /* PRECICE_IO_EXPORT_HPP_ */
