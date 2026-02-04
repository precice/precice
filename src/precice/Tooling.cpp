#include "precice/Tooling.hpp"

#include "fmt/color.h"
#include "precice/config/Configuration.hpp"
#include "precice/impl/versions.hpp"
#include "utils/Parallel.hpp"
#include "utils/Petsc.hpp"
#include "utils/fmt.hpp"
#include "xml/Printer.hpp"

namespace precice {

namespace tooling {

void printConfigReference(std::ostream &out, ConfigReferenceType reftype)
{
  precice::config::Configuration config;
  switch (reftype) {
  case ConfigReferenceType::XML:
    precice::xml::toDocumentation(out, config.getXMLTag());
    return;
  case ConfigReferenceType::DTD:
    precice::xml::toDTD(out, config.getXMLTag());
    return;
  case ConfigReferenceType::MD:
    out << "<!-- generated with preCICE " PRECICE_VERSION " -->\n";
    precice::xml::toMarkdown(out, config.getXMLTag());
    return;
  }
}

void checkConfiguration(const std::string &filename, const std::string &participant, int size)
{
  logging::setupLogging({[] {
    logging::BackendConfiguration config;

    // Console output
    config.format = "%ColorizedSeverity%%Message%";
    config.filter = "%Severity% >= info";
    config.type   = "stream";
    config.output = "stdout";

    return config;
  }()});

  // Lock the logging configuration to prevent the parser from changing it
  logging::lockConf();

  fmt::print("Checking {} for syntax and basic setup issues...\n", filename);
  config::Configuration config;
  logging::setMPIRank(0);
  const auto wasInitialized = utils::Parallel::isMPIInitialized();
  if (!wasInitialized) {
    utils::Parallel::initializeTestingMPI(nullptr, nullptr);
  }
  xml::ConfigurationContext context{
      participant,
      0,
      size};
  xml::configure(config.getXMLTag(), context, filename);
  fmt::print(fmt::emphasis::bold | fg(fmt::color::green), "No major issues detected\n");
  if (!wasInitialized) {
    utils::Petsc::finalize();
    utils::Parallel::finalizeTestingMPI();
  }
}

} // namespace tooling

std::string getVersionInformation()
{
  return {precice::versionInformation};
}

} // namespace precice
