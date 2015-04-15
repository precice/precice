// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_COUPLINGSCHEMECONFIGURATION_HPP_
#define PRECICE_CPLSCHEME_COUPLINGSCHEMECONFIGURATION_HPP_

#include "cplscheme/CouplingScheme.hpp"
#include "cplscheme/MultiCouplingScheme.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "cplscheme/impl/SharedPointer.hpp"
#include "cplscheme/Constants.hpp"
#include "mesh/SharedPointer.hpp"
#include "m2n/M2N.hpp"
#include "m2n/config/M2NConfiguration.hpp"
#include "precice/config/SharedPointer.hpp"
#include "utils/xml/XMLTag.hpp"
#include "tarch/logging/Log.h"
#include <vector>
#include <string>
#include "boost/tuple/tuple.hpp"
#include "precice/impl/MeshContext.hpp"

namespace precice {
  namespace cplscheme {
    class CompositionalCouplingScheme;
    class BaseCouplingScheme;
    namespace tests {
      class SerialImplicitCouplingSchemeTest;
      class ParallelImplicitCouplingSchemeTest;
    }
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace cplscheme {

/**
 * @brief Configuration for coupling schemes.
 */
class CouplingSchemeConfiguration : public utils::XMLTag::Listener
{
public:

  /**
   * @brief Constructor.
   *
   * @param parent [IN] Used to add subtags to hierarchical XML structure.
   * @param meshConfig [IN] For checking if a used mesh is defined.
   * @param comConfig [IN] For checking if a communication between participants
   *        to be coupled is defined.
   */
  CouplingSchemeConfiguration (
    utils::XMLTag&                            parent,
    const mesh::PtrMeshConfiguration&         meshConfig,
    const m2n::M2NConfiguration::SharedPointer&           m2nConfig);

  /**
   * @brief Destructor, empty.
   */
  virtual ~CouplingSchemeConfiguration() {}

  /**
   * @brief Check, if a coupling scheme is configured for a participant.
   */
  bool hasCouplingScheme ( const std::string& participantName ) const;

  /**
   * @brief Returns the configured (to be checked by isValid()) coupling scheme.
   */
  const PtrCouplingScheme& getCouplingScheme ( const std::string& participantName ) const;

  /**
   * @brief Returns the name of one dataset exchanged in the coupling scheme.
   */
  const std::string& getDataToExchange ( int index ) const;

  /**
   * @brief Callback method required when using utils::XMLTag.
   */
  virtual void xmlTagCallback ( utils::XMLTag& callingTag );

  /**
   * @brief Callback method required when using utils::XMLTag.
   */
  virtual void xmlEndTagCallback ( utils::XMLTag& callingTag );

  /**
   * @brief Adds a manually configured coupling scheme for a participant.
   */
  void addCouplingScheme (
    PtrCouplingScheme  cplScheme,
    const std::string& participantName );


private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  const std::string TAG;
  const std::string TAG_PARTICIPANTS;
  const std::string TAG_PARTICIPANT;
  const std::string TAG_EXCHANGE;
  const std::string TAG_MAX_TIME;
  const std::string TAG_MAX_TIMESTEPS;
  const std::string TAG_TIMESTEP_LENGTH;
  const std::string TAG_ABS_CONV_MEASURE;
  const std::string TAG_REL_CONV_MEASURE;
  const std::string TAG_RES_REL_CONV_MEASURE;
  const std::string TAG_MIN_ITER_CONV_MEASURE;
  const std::string TAG_MAX_ITERATIONS;
  const std::string TAG_CHECKPOINT;
  const std::string TAG_EXTRAPOLATION;

  const std::string ATTR_DATA;
  const std::string ATTR_MESH;
  const std::string ATTR_PARTICIPANT;
  const std::string ATTR_INITIALIZE;
  const std::string ATTR_TYPE;
  const std::string ATTR_FIRST;
  const std::string ATTR_SECOND;
  const std::string ATTR_VALUE;
  const std::string ATTR_VALID_DIGITS;
  const std::string ATTR_METHOD;
  const std::string ATTR_LIMIT;
  const std::string ATTR_MIN_ITERATIONS;
  const std::string ATTR_NAME;
  const std::string ATTR_TIMESTEP_INTERVAL;
  const std::string ATTR_FROM;
  const std::string ATTR_TO;
  const std::string ATTR_SUFFICES;
  const std::string ATTR_CONTROL;

  const std::string VALUE_SERIAL_EXPLICIT;
  const std::string VALUE_PARALLEL_EXPLICIT;
  const std::string VALUE_SERIAL_IMPLICIT;
  const std::string VALUE_PARALLEL_IMPLICIT;
  const std::string VALUE_MULTI;
  const std::string VALUE_UNCOUPLED;
  const std::string VALUE_FIXED;
  const std::string VALUE_FIRST_PARTICIPANT;

  struct Config
  {
    std::string type;
    std::string name;
    int checkpointTimestepInterval;
    std::vector<std::string> participants;
    std::string controller;
    bool setController;
    double maxTime;
    int maxTimesteps;
    double timestepLength;
    int validDigits;
    constants::TimesteppingMethod dtMethod;
    // @brief Tuples of exchange data, mesh, and participant name.
    typedef boost::tuple<mesh::PtrData, mesh::PtrMesh,std::string, std::string,bool> Exchange;
    std::vector<Exchange> exchanges;
    // @brief Tuples of data ID, mesh ID, and convergence measure.
    std::vector<boost::tuple<int,bool,std::string,impl::PtrConvergenceMeasure> > convMeasures;
    int maxIterations;
    int extrapolationOrder;

    Config()
    :
      type ( "" ),
      name ( "" ),
      checkpointTimestepInterval ( -1 ),
      participants (),
      controller ( "" ),
      setController( false ),
      maxTime ( CouplingScheme::UNDEFINED_TIME ),
      maxTimesteps ( CouplingScheme::UNDEFINED_TIMESTEPS ),
      timestepLength ( CouplingScheme::UNDEFINED_TIMESTEP_LENGTH ),
      validDigits ( 16 ),
      dtMethod ( constants::FIXED_DT ),
      exchanges (),
      convMeasures (),
      maxIterations ( -1 ),
      extrapolationOrder ( 0 )
    {}

  } _config;

  mesh::PtrMeshConfiguration _meshConfig;

  m2n::M2NConfiguration::SharedPointer _m2nConfig;

  PtrPostProcessingConfiguration _postProcConfig;

  bool _isValid;

  // @brief Map from participant name to coupling scheme (composition).
  std::map<std::string,PtrCouplingScheme> _couplingSchemes;

  // @brief If a participant has more than one coupling scheme, a composition is created.
  std::map<std::string,CompositionalCouplingScheme*> _couplingSchemeCompositions;

  void addTypespecifcSubtags (
    const std::string& type,
    utils::XMLTag&     tag );

  void addTagCheckpoint ( utils::XMLTag& tag );

  void addTransientLimitTags ( utils::XMLTag& tag );

  void addTagParticipants ( utils::XMLTag& tag );

  void addTagParticipant ( utils::XMLTag& tag );

  void addTagExchange ( utils::XMLTag& tag );

  void addTagAbsoluteConvergenceMeasure ( utils::XMLTag& tag );

  void addTagRelativeConvergenceMeasure ( utils::XMLTag& tag );

  void addTagResidualRelativeConvergenceMeasure ( utils::XMLTag& tag );

  void addTagMinIterationConvergenceMeasure ( utils::XMLTag& tag );

  void addBaseAttributesTagConvergenceMeasure ( utils::XMLTag& tag );

  void addTagMaxIterations ( utils::XMLTag& tag );

  void addTagExtrapolation ( utils::XMLTag& tag );

  void addTagPostProcessing ( utils::XMLTag& tag );

  void addAbsoluteConvergenceMeasure (
    const std::string & dataName,
    const std::string & meshName,
    double              limit,
    bool                suffices );

  void addRelativeConvergenceMeasure (
    const std::string & dataName,
    const std::string & meshName,
    double              limit,
    bool                suffices );

  void addResidualRelativeConvergenceMeasure (
    const std::string & dataName,
    const std::string & meshName,
    double              limit,
    bool                suffices );

  void addMinIterationConvergenceMeasure (
    const std::string & dataName,
    const std::string & meshName,
    int                 minIterations,
    bool                suffices );

  mesh::PtrData getData (
    const std::string & dataName,
    const std::string & meshName ) const;

  PtrCouplingScheme createSerialExplicitCouplingScheme (
    const std::string & accessor ) const;

  PtrCouplingScheme createParallelExplicitCouplingScheme (
      const std::string & accessor ) const;

  PtrCouplingScheme createSerialImplicitCouplingScheme (
    const std::string & accessor ) const;

  PtrCouplingScheme createParallelImplicitCouplingScheme (
    const std::string & accessor ) const;

  PtrCouplingScheme createMultiCouplingScheme (
      const std::string & accessor ) const;

  /*
   * @brief returns name of the actual scheme holder (i.e. server name)
   */
  std::string determineCouplingSchemeHolder (
    const std::string & accessorName ) const;

  constants::TimesteppingMethod getTimesteppingMethod (
    const std::string& method ) const;

  /**
   * @brief Adds configured exchange data to be sent or received to scheme.
   */
  void addDataToBeExchanged(
    BaseCouplingScheme& scheme,
    const std::string&  accessor) const;

  /**
   * @brief Adds configured exchange data to be sent or received to scheme.
   * Only used specifically for MultiCouplingScheme
   */
  void addMultiDataToBeExchanged(
    MultiCouplingScheme& scheme,
    const std::string&  accessor) const;

  void checkIfDataIsExchanged(
    int dataID) const;


  friend class tests::SerialImplicitCouplingSchemeTest; // For whitebox tests
  friend class tests::ParallelImplicitCouplingSchemeTest; // For whitebox tests

};

}} // namespace precice, cplscheme

#endif /* COUPLINGSCHEMECONFIGURATION_HPP_ */
