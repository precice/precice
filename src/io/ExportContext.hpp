#ifndef PRECICE_IO_EXPORTCONTEXT_HPP_
#define PRECICE_IO_EXPORTCONTEXT_HPP_

#include <memory>
#include <string>
#include "io/Export.hpp"

namespace precice::io {

struct ExportContext {
  // @brief Exporters performing the actual export.
  std::unique_ptr<Export> exporter;

  // @brief Path to export location.
  std::string location;

  // @brief Name of the mesh to export.
  std::string meshName;

  // @brief Exporting every N time windows (equals -1 when not set).
  int everyNTimeWindows = -1;

  // @brief If true, export is done in every iteration (also implicit).
  bool everyIteration = false;

  // @brief If true, updates the series file after every export. Otherwise only at the end
  bool updateSeries = false;

  // @brief type of the exporter (e.g. vtk).
  std::string type;
};

} // namespace precice::io

#endif /* PRECICE_IO_EXPORTCONTEXT_HPP_ */
