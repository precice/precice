#include "profiling/config/ProfilingConfiguration.hpp"
#include <cstdlib>
#include <filesystem>
#include <string_view>
#include "logging/LogMacros.hpp"
#include "profiling/EventUtils.hpp"
#include "utils/assertion.hpp"
#include "xml/ConfigParser.hpp"
#include "xml/XMLAttribute.hpp"
#include "xml/XMLTag.hpp"

bool precice::syncMode = false;

namespace precice::profiling {

namespace {
profiling::Mode fromString(std::string_view mode)
{
  if (mode == MODE_OFF) {
    return profiling::Mode::Off;
  } else if (mode == MODE_FUNDAMENTAL) {
    return profiling::Mode::Fundamental;
  } else if (mode == MODE_ALL) {
    return profiling::Mode::All;
  } else {
    PRECICE_UNREACHABLE("Unknown mode \"{}\"", mode);
  }
}
} // namespace

ProfilingConfiguration::ProfilingConfiguration(xml::XMLTag &parent)
{
  using namespace xml;

  XMLTag tag(*this, "profiling", XMLTag::OCCUR_NOT_OR_ONCE);
  tag.setDocumentation("Allows configuring the profiling functionality of preCICE.");

  auto attrMode = makeXMLAttribute<std::string>("mode", DEFAULT_MODE)
                      .setOptions({MODE_ALL, MODE_FUNDAMENTAL, MODE_OFF})
                      .setDocumentation("Operational modes of the profiling. "
                                        "\"fundamental\" will only write fundamental events. "
                                        "\"all\" writes all events.");
  tag.addAttribute(attrMode);

  auto attrFlush = makeXMLAttribute<int>("flush-every", DEFAULT_SYNC_EVERY)
                       .setDocumentation("Set the amount of event records that should be kept in memory before flushing them to file. "
                                         "One event consists out of multiple records."
                                         "0 keeps all records in memory and writes them at the end of the program, useful for slower network filesystems. "
                                         "1 writes records directly to the file, useful to get profiling data despite program crashes. "
                                         "Settings greater than 1 keep records in memory and write them to file in blocks, which is recommended.");
  tag.addAttribute(attrFlush);

  auto attrDirectory = makeXMLAttribute<std::string>("directory", DEFAULT_DIRECTORY)
                           .setDocumentation("Directory to use as a root directory to  write the events to. "
                                             "Events will be written to `<directory>/precice-profiling/`");
  tag.addAttribute(attrDirectory);

  auto attrSynchronize = xml::makeXMLAttribute("synchronize", false)
                             .setDocumentation("Enables additional inter- and intra-participant synchronization points. "
                                               "This avoids measuring blocking time for communication and other collective operations.");
  tag.addAttribute(attrSynchronize);

  parent.addSubtag(tag);
}

void ProfilingConfiguration::xmlTagCallback(
    const xml::ConfigurationContext &context,
    xml::XMLTag &                    tag)
{
  precice::syncMode = tag.getBooleanAttributeValue("synchronize");
  auto mode         = tag.getStringAttributeValue("mode");
  auto flushEvery   = tag.getIntAttributeValue("flush-every");
  auto directory    = std::filesystem::path(tag.getStringAttributeValue("directory"));
  PRECICE_CHECK(flushEvery >= 0, "You configured the profiling to flush-every=\"{}\", which is invalid. "
                                 "Please choose a number >= 0.");

  using namespace precice;
  auto &er = profiling::EventRegistry::instance();

  er.setWriteQueueMax(flushEvery);

  directory /= "precice-profiling";
  er.setDirectory(directory.string());

  er.setMode(fromString(mode));
}

void applyDefaults()
{
  precice::syncMode = false;
  auto &er          = profiling::EventRegistry::instance();

  er.setWriteQueueMax(DEFAULT_SYNC_EVERY);

  auto directory = std::filesystem::path(DEFAULT_DIRECTORY);
  directory /= "precice-profiling";
  er.setDirectory(directory.string());

  er.setMode(fromString(DEFAULT_MODE));
}

} // namespace precice::profiling
