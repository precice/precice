#pragma once

#include <string>

namespace precice {
namespace mesh {
class Mesh;
}
} // namespace precice

namespace precice {
namespace io {

/// Abstract base class of all classes exporting container data structures.
class Export {
public:
  Export &operator=(Export &&) = delete;

  virtual ~Export() {}

  /// Returns the export type ID.
  virtual int getType() const = 0;

  /**
   * @brief Does export. Has to be implemented in subclass.
   *
   * @param[in] name Filename (without path).
   * @param[in] location Location (path without filename).
   * @param[in] mesh Mesh to be exported.
   */
  virtual void doExport(
      const std::string &name,
      const std::string &location,
      mesh::Mesh &       mesh) = 0;
};

} // namespace io
} // namespace precice
