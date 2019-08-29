#pragma once

#include <memory>

namespace precice {
namespace io {

class Export;
class ExportConfiguration;

using PtrExport              = std::shared_ptr<Export>;
using PtrExportConfiguration = std::shared_ptr<ExportConfiguration>;

} // namespace io
} // namespace precice
