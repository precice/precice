#ifndef PRECICE_IO_EXPORTCONFIGURATION_HPP_
#define PRECICE_IO_EXPORTCONFIGURATION_HPP_

#include "io/ExportContext.hpp"
#include "io/Constants.hpp"
#include "io/SharedPointer.hpp"
#include "xml/XMLTag.hpp"
#include "logging/Logger.hpp"
#include "precice/Constants.hpp"
#include <string>
#include <list>

namespace precice {
namespace io {

/**
 * @brief Configuration class for exports.
 */
class ExportConfiguration : public utils::XMLTag::Listener
{
public:

  /**
   * @brief Name of the xml-tag corresponding to the ExportConfiguration.
   */
  //static const std::string& getTag();

  /**
   * @brief Constructor.
   */
  ExportConfiguration ( utils::XMLTag& parent );

  /**
   * @brief Parses the export configuration xml-tag.
   *
   * Requirements:
   * - xmlReader has to point to the tag corresponding to ExportConfiguration
   */
  //bool parseSubtag ( utils::XMLTag::XMLReader* xmlReader );

  /**
   * @brief Returns true, if xml-tag has been parsed successfully.
   */
  //bool isValid() const { return _isValid; }

  /**
   * @brief Returns the configured export context, valid if isValid() is true.
   */
  std::list<ExportContext>& exportContexts() { return _contexts; }

  /**
   * @brief Callback function required for use of automatic configuration.
   *
   * Is called by utils::XMLTag on automatic configuration every time an xml
   * tag and its attributes have been read.
   * @param callingTag [IN] XML tag currently read.
   * @param xmlReader  [IN] XML Reader responsible for reading the tag.
   * @return True, if the corresponding actions could be successfully performed.
   */
  virtual void xmlTagCallback ( utils::XMLTag& callingTag );

  /**
   * @brief Callback from automatic configuration. Not utilitzed here.
   */
  virtual void xmlEndTagCallback ( utils::XMLTag& callingTag ) {}

  void resetExports() { _contexts.clear(); }

private:

  // @brief Logging device.
  static logging::Logger _log;

  const std::string TAG;

  const std::string ATTR_LOCATION;
  const std::string ATTR_TYPE;
  const std::string ATTR_AUTO;
  const std::string VALUE_VTK;
  const std::string VALUE_VRML;

  const std::string ATTR_TIMESTEP_INTERVAL;
  const std::string ATTR_NEIGHBORS;
  const std::string ATTR_TRIGGER_SOLVER;
  const std::string ATTR_NORMALS;
  const std::string ATTR_EVERY_ITERATION;

  // @brief Flag indicating success of configuration.
  //bool _isValid;

  std::list<ExportContext> _contexts;
};

}} // namespace precice, io

#endif /* PRECICE_IO_EXPORTCONFIGURATION_HPP_ */
