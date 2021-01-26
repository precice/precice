#ifndef PRECICE_IO_EXPORTCONTEXT_HPP_
#define PRECICE_IO_EXPORTCONTEXT_HPP_

#include <string>
#include "io/SharedPointer.hpp"

namespace precice {
namespace io {

struct ExportContext {
  // @brief Exporters performing the actual export.
  io::PtrExport exporter;

  // @brief Path to export location.
  std::string location;

  // @brief Exporting every N time windows (equals -1 when not set).
  int everyNTimeWindows;

  // @brief If true, export is done in every iteration (also implicit).
  bool everyIteration;

  // @brief type of the exporter (e.g. vtk).
  std::string type;

  // @brief If true, normals are plotted.
  bool plotNormals;

  /**
   * @brief Constructor.
   */
  ExportContext()
      : exporter(),
        location(),
        everyNTimeWindows(-1),
        everyIteration(false),
        type(),
        plotNormals(false) {}
};

} // namespace io
} // namespace precice

#endif /* PRECICE_IO_EXPORTCONTEXT_HPP_ */
