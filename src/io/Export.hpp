#pragma once

#include <string>
#include <string_view>

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

  enum struct ExportKind : bool {
    TimeWindows,
    Iterations
  };

  Export(
      std::string_view  participantName,
      std::string_view  location,
      const mesh::Mesh &mesh,
      ExportKind        kind,
      int               frequency,
      int               rank,
      int               size)
      : _participantName(participantName),
        _location(location),
        _mesh(&mesh),
        _kind(kind),
        _frequency(frequency),
        _rank(rank),
        _size(size){};

  Export(const Export &) = delete;
  Export(Export &&)      = delete;
  Export &operator=(const Export &) = delete;
  Export &operator=(Export &&) = delete;

  /**
   * @brief Export the mesh and writes files.
   *
   * @param[in] index the index of this iteration or time window
   * @param[in] time the associated time of this time window
   */
  virtual void doExport(int index, double time) = 0;

protected:
  bool isParallel() const
  {
    return _size > 1;
  };
  std::string_view kindPrefix() const
  {
    return (_kind == ExportKind::TimeWindows) ? "dt" : "it";
  };

  std::string             _participantName;
  std::string             _location;
  const mesh::Mesh *const _mesh;
  ExportKind              _kind;
  int                     _frequency;
  int                     _rank;
  int                     _size;
};

} // namespace io
} // namespace precice
