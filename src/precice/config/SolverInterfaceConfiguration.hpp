#pragma once

#include <string>
#include "SharedPointer.hpp"
#include "boost/smart_ptr.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "m2n/M2N.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "mapping/SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"
#include "precice/impl/Participant.hpp"
#include "precice/impl/SharedPointer.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace config {

/**
 * @brief Configures class SolverInterfaceImpl from XML.
 */
class SolverInterfaceConfiguration : public xml::XMLTag::Listener {
public:
  SolverInterfaceConfiguration(xml::XMLTag &parent);

  /**
   * @brief Destructor.
   *
   * Deletes geometry configuration and coupling scheme configuration.
   */
  virtual ~SolverInterfaceConfiguration() {}

  /**
   * @brief Callback method required when using xml::XMLTag.
   */
  virtual void xmlTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &callingTag);

  /**
   * @brief Callback method required when using xml::XMLTag.
   */
  virtual void xmlEndTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &callingTag);

  /**
   * @brief Returns number of spatial dimensions configured.
   */
  int getDimensions() const;

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

  const PtrParticipantConfiguration &getParticipantConfiguration() const;

  const cplscheme::PtrCouplingSchemeConfiguration
  getCouplingSchemeConfiguration() const
  {
    return _couplingSchemeConfiguration;
  }

  /**
   * @brief For manual configuration in test cases.
   */
  void setDataConfiguration(mesh::PtrDataConfiguration config)
  {
    _dataConfiguration = config;
  }

  /**
   * @brief For manual configuration in test cases.
   */
  void setMeshConfiguration(mesh::PtrMeshConfiguration config)
  {
    _meshConfiguration = config;
  }

  /**
    * @brief For manual configuration in test cases.
    */
  void setParticipantConfiguration(PtrParticipantConfiguration config)
  {
    _participantConfiguration = config;
  }

private:
  logging::Logger _log{"config::SolverInterfaceConfiguration"};

  /// Spatial dimension of problem to be solved. Either 2 or 3.
  int _dimensions = -1;

  // @brief Participating solvers in the coupled simulation.
  //std::vector<impl::PtrParticipant> _participants;

  // @brief Index (in _participants) of solver accessing the interface.
  //int _indexAccessor;

  mesh::PtrDataConfiguration _dataConfiguration;

  mesh::PtrMeshConfiguration _meshConfiguration;

  m2n::M2NConfiguration::SharedPointer _m2nConfiguration;

  PtrParticipantConfiguration _participantConfiguration;

  cplscheme::PtrCouplingSchemeConfiguration _couplingSchemeConfiguration;
};

} // namespace config
} // namespace precice
