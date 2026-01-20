#pragma once

#include <string>
#include "logging/Logger.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
/// Enabled further inter- and intra-solver synchronisation
extern bool syncMode;
} // namespace precice

namespace precice::profiling {

constexpr int         DEFAULT_SYNC_EVERY = 50;
constexpr const char *DEFAULT_MODE       = "fundamental";
constexpr const char *DEFAULT_DIRECTORY  = ".";
constexpr const char *MODE_OFF           = "off";
constexpr const char *MODE_FUNDAMENTAL   = "fundamental";
constexpr const char *MODE_API           = "api";
constexpr const char *MODE_ALL           = "all";

/**
 * @brief Configuration class for exports.
 */
class ProfilingConfiguration final : public xml::XMLTag::Listener {
public:
  ProfilingConfiguration(xml::XMLTag &parent);

  ~ProfilingConfiguration() override = default;

  void xmlTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &callingTag) override;

  void xmlEndTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &callingTag) override {};

private:
  logging::Logger _log{"profiling::ProfilingConfiguration"};
};

void applyDefaults();

} // namespace precice::profiling
