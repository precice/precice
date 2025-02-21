#pragma once

#include <string>
#include <string_view>
#include <vector>

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

  virtual void exportSeries() const = 0;

protected:
  bool isParallel() const;

  std::string formatIndex(int index) const;

  bool keepExport(int index) const
  {
    return (_kind == ExportKind::Iterations) || ((_frequency > 0) && (index % _frequency == 0));
  }

  std::string             _participantName;
  std::string             _location;
  const mesh::Mesh *const _mesh;
  ExportKind              _kind;
  int                     _frequency;
  int                     _rank;
  int                     _size;

  struct Record {
    std::string filename;
    double      time;
  };

  std::vector<Record> _records;

  void writeSeriesFile(std::string_view filename) const;

  void recordExport(std::string filename, double time);
};

} // namespace io
} // namespace precice
