#pragma once

#include <memory>

namespace precice::io {

class Export;
class ExportConfiguration;

using PtrExport              = std::shared_ptr<Export>;
using PtrExportConfiguration = std::shared_ptr<ExportConfiguration>;

} // namespace precice::io
