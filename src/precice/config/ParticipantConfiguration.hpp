#ifndef PRECICE_CONFIG_PARTICIPANTCONFIGURATION_HPP_
#define PRECICE_CONFIG_PARTICIPANTCONFIGURATION_HPP_

#include "precice/impl/Participant.hpp"
#include "mesh/SharedPointer.hpp"
#include "mapping/SharedPointer.hpp"
#include "io/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "utils/xml/XMLTag.hpp"
#include <string>

namespace precice {
namespace config {

/**
 * @brief Performs XML configuration of a participant.
 */
class ParticipantConfiguration : public utils::XMLTag::Listener
{
public:

  // @brief Name of xml tag for this class in configuration file
  //static const std::string& getTag();

  /**
   * @brief Constructor.
   */
  ParticipantConfiguration (
    utils::XMLTag&                              parent,
    const mesh::PtrMeshConfiguration&           meshConfiguration);

  void setDimensions ( int dimensions );

  /**
   * @brief Reads the information parsed from an xml-file.
   */
  //bool parseSubtag ( utils::XMLTag::XMLReader* xmlReader );

  /**
   * @brief Callback function required for use of automatic configuration.
   *
   * @return True, if successful.
   */
  virtual void xmlTagCallback ( utils::XMLTag& callingTag );

  /**
   * @brief Callback function required for use of automatic configuration.
   *
   * @return True, if successful.
   */
  virtual void xmlEndTagCallback ( utils::XMLTag& callingTag );

  /**
   * @returns Returns true, if the xml-file parsing was successful.
   */
  //bool isValid() const;

  /**
   * @brief For manual configuration.
   */
  void addParticipant (
    const impl::PtrParticipant&             participant,
    const mapping::PtrMappingConfiguration& mappingConfig );

  /// Returns all configured participants.
  const std::vector<impl::PtrParticipant>& getParticipants() const;

private:

  struct WatchPointConfig
  {
    std::string name;
    std::string nameMesh;
    Eigen::VectorXd coordinates;
  };

  static logging::Logger _log;

  const std::string TAG;
  const std::string TAG_WRITE;
  const std::string TAG_READ;
  const std::string TAG_DATA_ACTION;
  const std::string TAG_USE_MESH;
  const std::string TAG_WATCH_POINT;
  const std::string TAG_SERVER;
  const std::string TAG_MASTER;

  const std::string ATTR_NAME;
  const std::string ATTR_SOURCE_DATA;
  const std::string ATTR_TARGET_DATA;
  const std::string ATTR_TIMING;
  const std::string ATTR_LOCAL_OFFSET;
  const std::string ATTR_ACTION_TYPE;
  const std::string ATTR_FROM;
  const std::string ATTR_SAFETY_FACTOR;
  const std::string ATTR_DECOMPOSITION;
  const std::string ATTR_PROVIDE;
  const std::string ATTR_MESH;
  const std::string ATTR_COORDINATE;
  const std::string ATTR_COMMUNICATION;
  const std::string ATTR_CONTEXT;
  const std::string ATTR_NETWORK;
  const std::string ATTR_EXCHANGE_DIRECTORY;

  const std::string VALUE_PRE_FILTER_POST_FILTER;
  const std::string VALUE_BROADCAST_FILTER;

  const std::string VALUE_VTK;
  const std::string VALUE_VRML;

  int _dimensions;

  mesh::PtrMeshConfiguration _meshConfig;

  mapping::PtrMappingConfiguration _mappingConfig;

  action::PtrActionConfiguration _actionConfig;

  io::PtrExportConfiguration _exportConfig;

  //bool _isValid;

  std::vector<impl::PtrParticipant> _participants;

  std::vector<WatchPointConfig> _watchPointConfigs;

  partition::ReceivedPartition::GeometricFilter getGeoFilter(const std::string& geoFilter) const;

  mesh::PtrMesh copy ( const mesh::PtrMesh& mesh ) const;

  const mesh::PtrData& getData (
    const mesh::PtrMesh& mesh,
    const std::string&   nameData ) const;

  mapping::PtrMapping getMapping ( const std::string& mappingName );

  void finishParticipantConfiguration ( const impl::PtrParticipant& participant );

};

}} // namespace precice, config

#endif /* PRECICE_CONFIG_PARTICIPANTCONFIGURATION_HPP_ */
