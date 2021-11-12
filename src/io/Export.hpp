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
  virtual ~Export() = default;

  Export()               = default;
  Export(const Export &) = delete;
  Export(Export &&)      = delete;
  Export &operator=(const Export &) = delete;
  Export &operator=(Export &&) = delete;

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
