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
#include "m2n/M2N.hpp"
#include "m2n/config/M2NConfiguration.hpp"
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

  /**
   * @brief Constructor.
   */
  SolverInterfaceConfiguration(utils::XMLTag& parent);

  /**
   * @brief Destructor.
   *
   * Deletes geometry configuration and coupling scheme configuration.
   */
  virtual ~SolverInterfaceConfiguration() {}

  /**
   * @brief Callback method required when using utils::XMLTag.
   */
  virtual void xmlTagCallback ( utils::XMLTag& callingTag );

  /**
   * @brief Callback method required when using utils::XMLTag.
   */
  virtual void xmlEndTagCallback ( utils::XMLTag& callingTag );

  /**
   * @brief Returns number of spatial dimensions configured.
   */
  int getDimensions() const;

  /**
   * @brief Returns true if preCICE is to be used as pure geometry interface.
   */
  bool isGeometryMode() const
  {
    return _geometryMode;
  }

  /**
   * @brief For manual configuration in test cases.
   */
  void setGeometryMode()
  {
    _geometryMode = true;
  }

  /**
   * @brief Returns true if a restart checkpoint should be read.
   */
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

  const m2n::M2NConfiguration::SharedPointer getM2NConfiguration() const
  {
    return _m2nConfiguration;
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

  /**
   * @brief For manual configuration in test cases.
   */
  void setDataConfiguration ( mesh::PtrDataConfiguration config )
  {
    _dataConfiguration = config;
  }

  /**
   * @brief For manual configuration in test cases.
   */
  void setMeshConfiguration ( mesh::PtrMeshConfiguration config )
  {
    _meshConfiguration = config;
  }

  /**
   * @brief For manual configuration in test cases.
   */
   void setGeometryConfiguration ( geometry::PtrGeometryConfiguration config )
   {
     _geometryConfiguration = config;
   }

   /**
    * @brief For manual configuration in test cases.
    */
   void setParticipantConfiguration ( PtrParticipantConfiguration config )
   {
     _participantConfiguration = config;
   }

private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  // Tag and subtag names used within this configuration.
  const std::string TAG;
  const std::string ATTR_DIMENSIONS;
  const std::string ATTR_GEOMETRY_MODE;
  const std::string ATTR_RESTART_MODE;
  const std::string ATTR_SPACETREE_NAME;

  // @brief Spatial dimension of problem to be solved. Either 2 or 3.
  int _dimensions;

  // @brief True, if pure geometry mode configuration without coupling is read.
  bool _geometryMode;

  // @brief True, if the coupled simulation is started from a checkpoint.
  bool _restartMode;

  // @brief Participating solvers in the coupled simulation.
  //std::vector<impl::PtrParticipant> _participants;

  // @brief Index (in _participants) of solver accessing the interface.
  //int _indexAccessor;

  mesh::PtrDataConfiguration _dataConfiguration;

  mesh::PtrMeshConfiguration _meshConfiguration;

  m2n::M2NConfiguration::SharedPointer _m2nConfiguration;

  geometry::PtrGeometryConfiguration _geometryConfiguration;

  spacetree::PtrSpacetreeConfiguration _spacetreeConfiguration;

  PtrParticipantConfiguration _participantConfiguration;

  cplscheme::PtrCouplingSchemeConfiguration _couplingSchemeConfiguration;
};

}} // namespace precice, config

#endif /* PRECICE_CONFIG_COUPLINGCONFIGURATION_HPP_ */
