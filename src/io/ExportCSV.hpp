#pragma once

#include <string>
#include "io/Export.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace io {

class ExportCSV : public Export {
public:
  ExportCSV(
      std::string_view  participantName,
      std::string_view  location,
      const mesh::Mesh &mesh,
      ExportKind        kind,
      int               frequency,
      int               rank,
      int               size);

  void doExport(int index, double time) final override;

  void exportSeries() const final override;

private:
  mutable logging::Logger _log{"io::ExportCSV"};
};

} // namespace io
} // namespace precice
