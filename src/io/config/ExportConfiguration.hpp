#pragma once

#include <list>
#include <string>
#include "io/ExportContext.hpp"
#include "logging/Logger.hpp"
#include "xml/XMLTag.hpp"

namespace precice::io {

/**
 * @brief Configuration class for exports.
 */
class ExportConfiguration : public xml::XMLTag::Listener {
public:
  ExportConfiguration(xml::XMLTag &parent);

  /**
   * @brief Returns the configured export context
   */
  std::list<ExportContext> &exportContexts()
  {
    return _contexts;
  }

  void xmlTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &callingTag) override;

  /// Callback from automatic configuration. Not utilitzed here.
  void xmlEndTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &callingTag) override {}

  void resetExports()
  {
    _contexts.clear();
  }

private:
  logging::Logger _log{"io::ExportConfiguration"};

  const std::string TAG = "export";

  const std::string ATTR_LOCATION = "directory";
  const std::string ATTR_TYPE     = "type";
  const std::string ATTR_AUTO     = "auto";
  const std::string VALUE_VTK     = "vtk";
  const std::string VALUE_VTU     = "vtu";
  const std::string VALUE_VTP     = "vtp";
  const std::string VALUE_CSV     = "csv";

  const std::string ATTR_EVERY_N_TIME_WINDOWS = "every-n-time-windows";
  const std::string ATTR_NEIGHBORS            = "neighbors";
  const std::string ATTR_EVERY_ITERATION      = "every-iteration";

  std::list<ExportContext> _contexts;
};

} // namespace precice::io
