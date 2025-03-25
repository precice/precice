#pragma once

#include <Eigen/Core>
#include <string>
#include <vector>
#include "action/SharedPointer.hpp"
#include "io/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "mapping/SharedPointer.hpp"
#include "mapping/config/MappingConfigurationTypes.hpp"
#include "mesh/SharedPointer.hpp"
#include "partition/ReceivedPartition.hpp"
#include "precice/impl/ParticipantState.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "time/Time.hpp"
#include "utils/networking.hpp"
#include "xml/XMLTag.hpp"

namespace precice::config {

/**
 * @brief Performs XML configuration of a participant.
 */
class ParticipantConfiguration : public xml::XMLTag::Listener {
public:
  ParticipantConfiguration(
      xml::XMLTag               &parent,
      mesh::PtrMeshConfiguration meshConfiguration);

  void setExperimental(bool experimental);
  void setRemeshing(bool allowed);

  /**
   * @brief Callback function required for use of automatic configuration.
   *
   * @return True, if successful.
   */
  void xmlTagCallback(
      const xml::ConfigurationContext &context,
      xml::XMLTag                     &callingTag) override;

  /**
   * @brief Callback function required for use of automatic configuration.
   *
   * @return True, if successful.
   */
  void xmlEndTagCallback(
      const xml::ConfigurationContext &context,
      xml::XMLTag                     &callingTag) override;

  /// Returns all configured participants.
  const std::vector<impl::PtrParticipant> &getParticipants() const;

  /// Returns a participant with the given name
  const impl::PtrParticipant getParticipant(const std::string &participantName) const;

  std::set<std::string> knownParticipants() const;

  bool hasParticipant(std::string_view name) const;

  std::string hintFor(std::string_view wrongName) const;

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
  const std::string TAG_PROVIDE_MESH   = "provide-mesh";
  const std::string TAG_RECEIVE_MESH   = "receive-mesh";
  const std::string TAG_WATCH_INTEGRAL = "watch-integral";
  const std::string TAG_WATCH_POINT    = "watch-point";
  const std::string TAG_INTRA_COMM     = "intra-comm";

  const std::string ATTR_NAME               = "name";
  const std::string ATTR_SOURCE_DATA        = "source-data";
  const std::string ATTR_TARGET_DATA        = "target-data";
  const std::string ATTR_TIMING             = "timing";
  const std::string ATTR_LOCAL_OFFSET       = "offset";
  const std::string ATTR_ACTION_TYPE        = "type";
  const std::string ATTR_FROM               = "from";
  const std::string ATTR_SAFETY_FACTOR      = "safety-factor";
  const std::string ATTR_GEOMETRIC_FILTER   = "geometric-filter";
  const std::string ATTR_DIRECT_ACCESS      = "direct-access";
  const std::string ATTR_API_ACCESS         = "api-access";
  const std::string ATTR_PROVIDE            = "provide";
  const std::string ATTR_MESH               = "mesh";
  const std::string ATTR_COORDINATE         = "coordinate";
  const std::string ATTR_COMMUNICATION      = "communication";
  const std::string ATTR_CONTEXT            = "context";
  const std::string ATTR_NETWORK            = "network";
  const std::string ATTR_EXCHANGE_DIRECTORY = "exchange-directory";
  const std::string ATTR_SCALE_WITH_CONN    = "scale-with-connectivity";

  const std::string VALUE_FILTER_ON_SECONDARY_RANKS = "on-secondary-ranks";
  const std::string VALUE_FILTER_ON_PRIMARY_RANK    = "on-primary-rank";
  const std::string VALUE_NO_FILTER                 = "no-filter";

  const std::string VALUE_VTK = "vtk";
  const std::string VALUE_VTU = "vtu";
  const std::string VALUE_VTP = "vtp";
  const std::string VALUE_CSV = "csv";

  bool _experimental = false;
  bool _remeshing    = false;

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
      const std::string   &nameData) const;

  mapping::PtrMapping getMapping(const std::string &mappingName);

  // Does this participant already define a primary tag?
  // This context information is needed in xmlEndTagCallback to create a default
  // primary com if required (i.e. no solution yet defined and parallel).
  bool _isIntraCommDefined = false;

  void finishParticipantConfiguration(
      const xml::ConfigurationContext &context,
      const impl::PtrParticipant      &participant);

  /// Check whether a mapping to the same mesh and with similar data fields already exists
  void checkIllDefinedMappings(
      const mapping::MappingConfiguration::ConfiguredMapping &mapping,
      const impl::PtrParticipant                             &participant);
};

} // namespace precice::config
