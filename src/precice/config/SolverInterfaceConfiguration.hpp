// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CONFIG_COUPLINGCONFIGURATION_HPP_
#define PRECICE_CONFIG_COUPLINGCONFIGURATION_HPP_

#include "SharedPointer.hpp"
#include "precice/impl/Participant.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"
#include "geometry/SharedPointer.hpp"
#include "spacetree/SharedPointer.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "mapping/SharedPointer.hpp"
#include "com/SharedPointer.hpp"
#include "tarch/logging/Log.h"
#include "utils/xml/XMLTag.hpp"
#include <string>
#include "boost/smart_ptr.hpp"

namespace precice {
namespace config {

/**
 * @brief Configures class SolverInterfaceImpl from XML.
 */
class SolverInterfaceConfiguration : public utils::XMLTag::Listener
{
public:

  // @brief Name of xml tag for this class in configuration file
  //static const std::string& getTag();

  /**
   * @brief Constructor.
   */
  SolverInterfaceConfiguration(utils::XMLTag& parent);

  /**
   * @brief Destructor.
   *
   * Deletes geometry configuration and coupling scheme configuration.
   */
  ~SolverInterfaceConfiguration() {}

  /**
   * @brief Reads the information parsed from an xml-file.
   */
  //bool parseSubtag ( utils::XMLTag::XMLReader* xmlReader );

  /**
   * @returns Returns true, if the xml-file parsing was successful.
   */
  //bool isValid() const;

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
   * @brief Adds UncoupledCouplingScheme in case of geometry mode.
   */
  virtual void xmlEndTagCallback ( utils::XMLTag& callingTag );

  /**
   * @brief Returns number of spatial dimensions configured.
   */
  int getDimensions() const;

  bool isGeometryMode() const
  {
    return _geometryMode;
  }

  void setGeometryMode()
  {
    _geometryMode = true;
  }

  bool isRestartMode() const
  {
    return _restartMode;
  }

  const mesh::PtrDataConfiguration getDataConfiguration() const
  {
    return _dataConfiguration;
  }

  const mesh::PtrMeshConfiguration getMeshConfiguration() const
  {
    return _meshConfiguration;
  }

  const com::PtrCommunicationConfiguration getCommunicationConfiguration() const
  {
    return _comConfiguration;
  }

  const geometry::PtrGeometryConfiguration getGeometryConfiguration() const
  {
    return _geometryConfiguration;
  }

  const spacetree::PtrSpacetreeConfiguration& getSpacetreeConfiguration() const;

  const PtrParticipantConfiguration& getParticipantConfiguration() const;

  const cplscheme::PtrCouplingSchemeConfiguration
  getCouplingSchemeConfiguration() const
  {
    return _couplingSchemeConfiguration;
  }

  void setDataConfiguration ( mesh::PtrDataConfiguration config )
  {
    _dataConfiguration = config;
  }

  void setMeshConfiguration ( mesh::PtrMeshConfiguration config )
  {
    _meshConfiguration = config;
  }

  /**
   * @brief For manual configuration.
   */
   void setGeometryConfiguration ( geometry::PtrGeometryConfiguration config )
   {
     _geometryConfiguration = config;
   }

   void setParticipantConfiguration ( PtrParticipantConfiguration config )
   {
     _participantConfiguration = config;
   }

   /**
    * @brief Is meant for test cases only!
    */
//   void setIsValid()
//   {
//     _isValid = true;
//   }

private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  const std::string TAG;
  const std::string ATTR_DIMENSIONS;
  const std::string ATTR_GEOMETRY_MODE;
  const std::string ATTR_RESTART_MODE;
  const std::string ATTR_SPACETREE_NAME;

  // @brief True, when a valid xml-tag has been parsed.
  //bool _isValid;

  // @brief Spatial dimension of problem to be solved. Either 2 or 3.
  int _dimensions;

  // @brief True, if pure geometry mode configuration without coupling is read.
  bool _geometryMode;

  // @brief True, if the coupled simulation is started from a checkpoint.
  bool _restartMode;

  // @brief Participating solvers in the coupled simulation.
  std::vector<impl::PtrParticipant> _participants;

  mesh::PtrDataConfiguration _dataConfiguration;

  mesh::PtrMeshConfiguration _meshConfiguration;

  com::PtrCommunicationConfiguration _comConfiguration;

  geometry::PtrGeometryConfiguration _geometryConfiguration;

  // @brief Spacetree name -> spacetree configuration for that geometry.
  spacetree::PtrSpacetreeConfiguration _spacetreeConfiguration;

  PtrParticipantConfiguration _participantConfiguration;

  // @brief Holds coupling information between participants.
  cplscheme::PtrCouplingSchemeConfiguration _couplingSchemeConfiguration;

  // @brief Index (in participant vector) of solver accessing the interface.
  int _indexAccessor;

//  int getDataDimensions ( const std::string& typeName );

//  int readDimensions ( utils::XMLTag::XMLReader* xmlReader );
//
//  bool readGeometryMode ( utils::XMLTag::XMLReader* xmlReader );
//
//  bool readRestartMode ( utils::XMLTag::XMLReader* xmlReader );
};

}} // namespace precice, config

#endif /* PRECICE_CONFIG_COUPLINGCONFIGURATION_HPP_ */
