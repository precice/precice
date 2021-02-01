#pragma once

#include <Eigen/Core>
#include <string>
#include <vector>
#include "action/SharedPointer.hpp"
#include "io/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "mapping/SharedPointer.hpp"
#include "mapping/config/MappingConfiguration.hpp"
#include "mesh/SharedPointer.hpp"
#include "partition/ReceivedPartition.hpp"
#include "precice/impl/Participant.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "utils/networking.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace config {

/**
 * @brief Performs XML configuration of a participant.
 */
class ParticipantConfiguration : public xml::XMLTag::Listener {
public:
  ParticipantConfiguration(
      xml::XMLTag &                     parent,
      const mesh::PtrMeshConfiguration &meshConfiguration);

  void setDimensions(int dimensions);

  /**
   * @brief Callback function required for use of automatic configuration.
   *
   * @return True, if successful.
   */
  virtual void xmlTagCallback(
      const xml::ConfigurationContext &context,
      xml::XMLTag &                    callingTag);

  /**
   * @brief Callback function required for use of automatic configuration.
   *
   * @return True, if successful.
   */
  virtual void xmlEndTagCallback(
      const xml::ConfigurationContext &context,
      xml::XMLTag &                    callingTag);

  /// Returns all configured participants.
  const std::vector<impl::PtrParticipant> &getParticipants() const;

private:
  struct WatchPointConfig {
    std::string     name;
    std::string     nameMesh;
    Eigen::VectorXd coordinates;
  };

  struct WatchIntegralConfig {
    std::string name;
    std::string nameMesh;
    bool        isScalingOn;
  };

  mutable logging::Logger _log{"config::ParticipantConfiguration"};

  const std::string TAG                = "participant";
  const std::string TAG_WRITE          = "write-data";
  const std::string TAG_READ           = "read-data";
  const std::string TAG_DATA_ACTION    = "data-action";
  const std::string TAG_USE_MESH       = "use-mesh";
  const std::string TAG_WATCH_INTEGRAL = "watch-integral";
  const std::string TAG_WATCH_POINT    = "watch-point";
  const std::string TAG_MASTER         = "master";

  const std::string ATTR_NAME               = "name";
  const std::string ATTR_SOURCE_DATA        = "source-data";
  const std::string ATTR_TARGET_DATA        = "target-data";
  const std::string ATTR_TIMING             = "timing";
  const std::string ATTR_LOCAL_OFFSET       = "offset";
  const std::string ATTR_ACTION_TYPE        = "type";
  const std::string ATTR_FROM               = "from";
  const std::string ATTR_SAFETY_FACTOR      = "safety-factor";
  const std::string ATTR_GEOMETRIC_FILTER   = "geometric-filter";
  const std::string ATTR_PROVIDE            = "provide";
  const std::string ATTR_MESH               = "mesh";
  const std::string ATTR_COORDINATE         = "coordinate";
  const std::string ATTR_COMMUNICATION      = "communication";
  const std::string ATTR_CONTEXT            = "context";
  const std::string ATTR_NETWORK            = "network";
  const std::string ATTR_EXCHANGE_DIRECTORY = "exchange-directory";
  const std::string ATTR_SCALE_WITH_CONN    = "scale-with-connectivity";

  const std::string VALUE_FILTER_ON_SLAVES = "on-slaves";
  const std::string VALUE_FILTER_ON_MASTER = "on-master";
  const std::string VALUE_NO_FILTER        = "no-filter";

  const std::string VALUE_VTK = "vtk";

  int _dimensions = 0;

  mesh::PtrMeshConfiguration _meshConfig;

  mapping::PtrMappingConfiguration _mappingConfig;

  action::PtrActionConfiguration _actionConfig;

  io::PtrExportConfiguration _exportConfig;

  std::vector<impl::PtrParticipant> _participants;

  std::vector<WatchPointConfig> _watchPointConfigs;

  std::vector<WatchIntegralConfig> _watchIntegralConfigs;

  partition::ReceivedPartition::GeometricFilter getGeoFilter(const std::string &geoFilter) const;

  mesh::PtrMesh copy(const mesh::PtrMesh &mesh) const;

  const mesh::PtrData &getData(
      const mesh::PtrMesh &mesh,
      const std::string &  nameData) const;

  mapping::PtrMapping getMapping(const std::string &mappingName);

  // Does this participant already define a master tag?
  // This context information is needed in xmlEndTagCallback to create a default
  // master com if required (i.e. no solution yet defined and parallel).
  bool _isMasterDefined = false;

  void finishParticipantConfiguration(
      const xml::ConfigurationContext &context,
      const impl::PtrParticipant &     participant);

  /// Check whether a mapping to the same mesh and with similar data fields already exists
  void checkIllDefinedMappings(
      const mapping::MappingConfiguration::ConfiguredMapping &mapping,
      const impl::PtrParticipant &                            participant);
};

} // namespace config
} // namespace precice
