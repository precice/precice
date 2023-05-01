#pragma once

#include <map>
#include <string>
#include <tuple>
#include <vector>

#include "acceleration/SharedPointer.hpp"
#include "cplscheme/Constants.hpp"
#include "cplscheme/CouplingScheme.hpp"
#include "cplscheme/MultiCouplingScheme.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "cplscheme/impl/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "m2n/config/M2NConfiguration.hpp"
// #include "mesh/GlobalData.hpp" // this probably shouldn't be needed, but without this compiler throws warnings.
#include "mesh/SharedPointer.hpp"
#include "precice/config/SharedPointer.hpp"
#include "precice/impl/MeshContext.hpp"
#include "precice/types.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace cplscheme {
class CompositionalCouplingScheme;
class BiCouplingScheme;
} // namespace cplscheme
} // namespace precice

// Forward declaration to friend the boost test struct
namespace CplSchemeTests {
namespace ParallelImplicitCouplingSchemeTests {
struct testParseConfigurationWithRelaxation;
}
namespace SerialImplicitCouplingSchemeTests {
struct testParseConfigurationWithRelaxation;
}
} // namespace CplSchemeTests

// ----------------------------------------------------------- CLASS DEFINITION
namespace precice {
namespace cplscheme {
class MultiCouplingScheme;

/// Configuration for coupling schemes.
class CouplingSchemeConfiguration : public xml::XMLTag::Listener {
public:
  /**
   * @brief Constructor.
   *
   * @param[in] parent  Used to add subtags to hierarchical XML structure.
   * @param[in] meshConfig For checking if a used mesh is defined.
   * @param[in] m2nConfig For checking if a communication between participants to be coupled is defined.
   * @param[in] participantConfig For checking waveform degree.
   */
  CouplingSchemeConfiguration(
      xml::XMLTag &                        parent,
      mesh::PtrMeshConfiguration           meshConfig,
      m2n::M2NConfiguration::SharedPointer m2nConfig,
      config::PtrParticipantConfiguration  participantConfig);

  void setExperimental(bool experimental);

  /// Destructor, empty.
  virtual ~CouplingSchemeConfiguration() {}

  /// Check, if a coupling scheme is configured for a participant.
  bool hasCouplingScheme(const std::string &participantName) const;

  /// Returns the configured coupling scheme.
  const PtrCouplingScheme &getCouplingScheme(const std::string &participantName) const;

  /// Returns the name of one dataset exchanged in the coupling scheme.
  const std::string &getDataToExchange(int index) const;

  /// Callback method required when using xml::XMLTag.
  virtual void xmlTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &callingTag);

  /// Callback method required when using xml::XMLTag.
  virtual void xmlEndTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &callingTag);

  /// Adds a manually configured coupling scheme for a participant.
  void addCouplingScheme(const PtrCouplingScheme &cplScheme, const std::string &participantName);

private:
  mutable logging::Logger _log{"cplscheme::CouplingSchemeConfiguration"};

  const std::string TAG;
  const std::string TAG_PARTICIPANTS;
  const std::string TAG_PARTICIPANT;
  const std::string TAG_EXCHANGE;
  const std::string TAG_MAX_TIME;
  const std::string TAG_MAX_TIME_WINDOWS;
  const std::string TAG_TIME_WINDOW_SIZE;
  const std::string TAG_ABS_CONV_MEASURE;
  const std::string TAG_REL_CONV_MEASURE;
  const std::string TAG_RES_REL_CONV_MEASURE;
  const std::string TAG_MIN_ITER_CONV_MEASURE;
  const std::string TAG_MAX_ITERATIONS;
  const std::string TAG_EXTRAPOLATION;

  const std::string ATTR_DATA;
  const std::string ATTR_MESH;
  const std::string ATTR_PARTICIPANT;
  const std::string ATTR_INITIALIZE;
  const std::string ATTR_EXCHANGE_SUBSTEPS;
  const std::string ATTR_TYPE;
  const std::string ATTR_FIRST;
  const std::string ATTR_SECOND;
  const std::string ATTR_VALUE;
  const std::string ATTR_VALID_DIGITS;
  const std::string ATTR_METHOD;
  const std::string ATTR_LIMIT;
  const std::string ATTR_MIN_ITERATIONS;
  const std::string ATTR_NAME;
  const std::string ATTR_FROM;
  const std::string ATTR_TO;
  const std::string ATTR_SUFFICES;
  const std::string ATTR_STRICT;
  const std::string ATTR_CONTROL;

  const std::string VALUE_SERIAL_EXPLICIT;
  const std::string VALUE_PARALLEL_EXPLICIT;
  const std::string VALUE_SERIAL_IMPLICIT;
  const std::string VALUE_PARALLEL_IMPLICIT;
  const std::string VALUE_MULTI;
  const std::string VALUE_FIXED;
  const std::string VALUE_FIRST_PARTICIPANT;

  struct ConvergenceMeasureDefintion {
    mesh::PtrData               data;
    bool                        suffices;
    bool                        strict;
    std::string                 meshName;
    impl::PtrConvergenceMeasure measure;
    bool                        doesLogging;
  };
  struct ConvergenceMeasureDefintionGlobalData {
    mesh::PtrData               globalData;
    bool                        suffices;
    bool                        strict;
    impl::PtrConvergenceMeasure measure;
    bool                        doesLogging;
  };

  struct Config {
    std::string                   type;
    std::string                   name;
    std::vector<std::string>      participants;
    std::string                   controller;
    bool                          setController  = false;
    double                        maxTime        = CouplingScheme::UNDEFINED_MAX_TIME;
    int                           maxTimeWindows = CouplingScheme::UNDEFINED_TIME_WINDOWS;
    double                        timeWindowSize = CouplingScheme::UNDEFINED_TIME_WINDOW_SIZE;
    int                           validDigits    = 16;
    constants::TimesteppingMethod dtMethod       = constants::FIXED_TIME_WINDOW_SIZE;

    struct Exchange {
      mesh::PtrData data;
      mesh::PtrMesh mesh;
      std::string   from;
      std::string   to;
      bool          requiresInitialization;
      bool          exchangeSubsteps;
    };
    struct GlobalExchange {
      mesh::PtrData globalData;
      std::string   from;
      std::string   to;
      bool          requiresInitialization;
    };
    std::vector<Exchange>                              exchanges;
    std::vector<GlobalExchange>                        globalExchanges;
    std::vector<ConvergenceMeasureDefintion>           convergenceMeasureDefinitions;
    std::vector<ConvergenceMeasureDefintionGlobalData> convergenceMeasureDefinitionsGlobalData;
    int                                                maxIterations      = -1;
    int                                                extrapolationOrder = 0;

    bool hasExchange(const Exchange &totest) const
    {
      return std::any_of(exchanges.begin(), exchanges.end(), [&totest](const auto &ex) {
        return ex.from == totest.from && ex.to == totest.to && ex.data->getName() == totest.data->getName() && ex.mesh->getName() == totest.mesh->getName();
      });
    }

    bool hasGlobalExchange(const GlobalExchange &totest) const
    {
      return std::any_of(globalExchanges.begin(), globalExchanges.end(), [&totest](const auto &gex) {
        return gex.from == totest.from && gex.to == totest.to && gex.globalData->getName() == totest.globalData->getName();
      });
    }
  } _config;

  mesh::PtrMeshConfiguration _meshConfig;

  m2n::M2NConfiguration::SharedPointer _m2nConfig;

  acceleration::PtrAccelerationConfiguration _accelerationConfig;

  precice::config::PtrParticipantConfiguration _participantConfig;

  /// Map from participant name to coupling scheme (composition).
  std::map<std::string, PtrCouplingScheme> _couplingSchemes;

  /// If a participant has more than one coupling scheme, a composition is created.
  std::map<std::string, CompositionalCouplingScheme *> _couplingSchemeCompositions;

  void addTypespecifcSubtags(const std::string &type, xml::XMLTag &tag);

  void addTransientLimitTags(const std::string &type, xml::XMLTag &tag);

  void addTagParticipants(xml::XMLTag &tag);

  void addTagParticipant(xml::XMLTag &tag);

  void addTagExchange(xml::XMLTag &tag);

  void addTagAbsoluteConvergenceMeasure(xml::XMLTag &tag);

  void addTagRelativeConvergenceMeasure(xml::XMLTag &tag);

  void addTagResidualRelativeConvergenceMeasure(xml::XMLTag &tag);

  void addTagMinIterationConvergenceMeasure(xml::XMLTag &tag);

  void addBaseAttributesTagConvergenceMeasure(xml::XMLTag &tag);

  void addTagMaxIterations(xml::XMLTag &tag);

  void addTagExtrapolation(xml::XMLTag &tag);

  void addTagAcceleration(xml::XMLTag &tag);

  void addAbsoluteConvergenceMeasure(
      const std::string &dataName,
      const std::string &meshName,
      double             limit,
      bool               suffices,
      bool               strict);

  void addAbsoluteConvergenceMeasureGlobalData(
      const std::string &dataName,
      double             limit,
      bool               suffices,
      bool               strict);

  void addRelativeConvergenceMeasure(
      const std::string &dataName,
      const std::string &meshName,
      double             limit,
      bool               suffices,
      bool               strict);

  void addRelativeConvergenceMeasureGlobalData(
      const std::string &dataName,
      double             limit,
      bool               suffices,
      bool               strict);

  void addResidualRelativeConvergenceMeasure(
      const std::string &dataName,
      const std::string &meshName,
      double             limit,
      bool               suffices,
      bool               strict);

  void addResidualRelativeConvergenceMeasureGlobalData(
      const std::string &dataName,
      double             limit,
      bool               suffices,
      bool               strict);

  void addMinIterationConvergenceMeasure(
      const std::string &dataName,
      const std::string &meshName,
      int                minIterations,
      bool               suffices,
      bool               strict);

  void addMinIterationConvergenceMeasureGlobalData(
      const std::string &dataName,
      int                minIterations,
      bool               suffices,
      bool               strict);

  mesh::PtrData getData(
      const std::string &dataName,
      const std::string &meshName) const;

  mesh::PtrData getGlobalData(
      const std::string &dataName) const;

  mesh::PtrData findDataByID(
      int ID) const;

  mesh::PtrData findGlobalDataByID(
      int ID) const;

  PtrCouplingScheme createSerialExplicitCouplingScheme(
      const std::string &accessor) const;

  PtrCouplingScheme createParallelExplicitCouplingScheme(
      const std::string &accessor) const;

  PtrCouplingScheme createSerialImplicitCouplingScheme(
      const std::string &accessor) const;

  PtrCouplingScheme createParallelImplicitCouplingScheme(
      const std::string &accessor) const;

  PtrCouplingScheme createMultiCouplingScheme(
      const std::string &accessor) const;

  constants::TimesteppingMethod getTimesteppingMethod(
      const std::string &method) const;

  /// Adds configured exchange data (mesh-associated as well as global) to be sent or received to scheme.
  void addDataToBeExchanged(
      BiCouplingScheme & scheme,
      const std::string &accessor) const;

  /**
   * @brief Adds configured exchange data to be sent or received to scheme.
   * Only used specifically for MultiCouplingScheme
   */
  void addMultiDataToBeExchanged(
      MultiCouplingScheme &scheme,
      const std::string &  accessor) const;

  void checkIfDataIsExchanged(
      DataID dataID) const;

  void checkIfGlobalDataIsExchanged(
      DataID dataID) const;

  void checkSerialImplicitAccelerationData(
      DataID dataID, const std::string &first, const std::string &second) const;

  void addConvergenceMeasures(
      BaseCouplingScheme *                            scheme,
      const std::string &                             participant,
      const std::vector<ConvergenceMeasureDefintion> &convergenceMeasureDefinitions) const;

  void addConvergenceMeasuresGlobalData(
      BaseCouplingScheme *                                      scheme,
      const std::string &                                       participant,
      const std::vector<ConvergenceMeasureDefintionGlobalData> &convergenceMeasureDefinitionsGlobalData) const;

  void setSerialAcceleration(
      BaseCouplingScheme *scheme,
      const std::string & first,
      const std::string & second) const;

  void setParallelAcceleration(
      BaseCouplingScheme *scheme,
      const std::string & participant) const;

  friend struct CplSchemeTests::ParallelImplicitCouplingSchemeTests::testParseConfigurationWithRelaxation; // For whitebox tests
  friend struct CplSchemeTests::SerialImplicitCouplingSchemeTests::testParseConfigurationWithRelaxation;   // For whitebox tests

  /**
   * @brief Helper function to check that waveform-degree and substep exchange are compatible.
   *
   * The following rules are checked:
   *
   * 1) If waveform-degree="0", then user must set substeps="false", because constant interpolation (zeroth degree) is intended for debugging and user should use first degree instead.
   * 2) If waveform-degree="1", then any configuration for substeps is allowed. The user might want to set substeps="false" for better performance.
   * 3) If waveform-degree="2" or greater, the user must set substeps="true", because subcycling and exchange of substeps is required for higher-degree B-splines.
   *
   * @param exchange The Exchange being checked.
   */
  void checkSubstepExchangeWaveformDegree(const Config::Exchange &exchange) const;
};
} // namespace cplscheme
} // namespace precice
