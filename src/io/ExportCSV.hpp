#pragma once

#include <string>
#include "io/Export.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace io {

class ExportCSV : public Export {
public:
  virtual void doExport(
      const std::string &name,
      const std::string &location,
      const mesh::Mesh & mesh);

private:
  mutable logging::Logger _log{"io::ExportCSV"};
};

} // namespace io
} // namespace precice
