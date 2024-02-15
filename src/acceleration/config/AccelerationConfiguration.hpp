#pragma once

#include <map>
#include <string>
#include <vector>
#include "acceleration/Acceleration.hpp"
#include "acceleration/IQNIMVJAcceleration.hpp"
#include "acceleration/SharedPointer.hpp"
#include "acceleration/impl/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"
#include "precice/config/SharedPointer.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace acceleration {

class AccelerationConfiguration : public xml::XMLTag::Listener {
public:
  AccelerationConfiguration(const mesh::PtrMeshConfiguration &meshConfig);

  /// Returns the configured coupling scheme.
  PtrAcceleration getAcceleration();

  /// Callback method required when using xml::XMLTag.
  virtual void xmlTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &callingTag);

  /// Callback method required when using xml::XMLTag.
  virtual void xmlEndTagCallback(const xml::ConfigurationContext &context, xml::XMLTag &callingTag);

  /// Removes configured acceleration.
  void clear();

  /// Connect tags.
  void connectTags(xml::XMLTag &tag);

  std::vector<std::string> &getNeededMeshes()
  {
    return _neededMeshes;
  }

private:
  logging::Logger _log{"acceleration::AccelerationConfiguration"};

  const std::string TAG;
  const std::string TAG_RELAX;
  const std::string TAG_INIT_RELAX;
  const std::string TAG_MAX_USED_ITERATIONS;
  const std::string TAG_TIME_WINDOWS_REUSED;
  const std::string TAG_DATA;
  const std::string TAG_FILTER;
  const std::string TAG_ESTIMATEJACOBIAN;
  const std::string TAG_PRECONDITIONER;
  const std::string TAG_IMVJRESTART;

  const std::string ATTR_NAME;
  const std::string ATTR_MESH;
  const std::string ATTR_SCALING;
  const std::string ATTR_VALUE;
  const std::string ATTR_MIN;
  const std::string ATTR_MAX;
  const std::string ATTR_ENFORCE;
  const std::string ATTR_SINGULARITYLIMIT;
  const std::string ATTR_TYPE;
  const std::string ATTR_BUILDJACOBIAN;
  const std::string ATTR_IMVJCHUNKSIZE;
  const std::string ATTR_RSLS_REUSED_TIME_WINDOWS;
  const std::string ATTR_RSSVD_TRUNCATIONEPS;
  const std::string ATTR_PRECOND_NONCONST_TIME_WINDOWS;

  const std::string VALUE_CONSTANT;
  const std::string VALUE_AITKEN;
  const std::string VALUE_IQNILS;
  const std::string VALUE_IQNIMVJ;
  const std::string VALUE_QR1FILTER;
  const std::string VALUE_QR1_ABSFILTER;
  const std::string VALUE_QR2FILTER;
  const std::string VALUE_CONSTANT_PRECONDITIONER;
  const std::string VALUE_VALUE_PRECONDITIONER;
  const std::string VALUE_RESIDUAL_PRECONDITIONER;
  const std::string VALUE_RESIDUAL_SUM_PRECONDITIONER;
  const std::string VALUE_LS_RESTART;
  const std::string VALUE_ZERO_RESTART;
  const std::string VALUE_SVD_RESTART;
  const std::string VALUE_SLIDE_RESTART;
  const std::string VALUE_NO_RESTART;

  const mesh::PtrMeshConfiguration _meshConfig;

  std::string _meshName;

  // acceleration method
  PtrAcceleration _acceleration;

  std::vector<std::string> _neededMeshes;

  impl::PtrPreconditioner _preconditioner;

  std::set<std::pair<std::string, std::string>> _uniqueDataAndMeshNames;

  struct ConfigurationData {
    std::vector<int>      dataIDs;
    std::map<int, double> scalings;
    std::map<int, double> lowerBounds;
    std::map<int, double> upperBounds;
    std::string           type;
    double                relaxationFactor           = 0;
    bool                  forceInitialRelaxation     = false;
    int                   maxIterationsUsed          = 0;
    int                   timeWindowsReused          = 0;
    int                   filter                     = Acceleration::NOFILTER;
    int                   imvjRestartType            = 0;
    int                   imvjChunkSize              = 0;
    int                   imvjRSLS_reusedTimeWindows = 0;
    int                   precond_nbNonConstTWindows = -1;
    double                singularityLimit           = 0;
    double                imvjRSSVD_truncationEps    = 0;
    bool                  estimateJacobian           = false;
    bool                  alwaysBuildJacobian        = false;
    std::string           preconditionerType;

    std::vector<double> scalingFactorsInOrder() const;
  } _config;

  const struct DefaultValuesIQN {
    double      relaxationFactor           = 0.1;
    int         maxIterationsUsed          = 100;
    int         timeWindowsReused          = 10;
    int         filter                     = Acceleration::QR2FILTER;
    double      singularityLimit           = 1e-2;
    std::string preconditionerType         = "residual-sum";
    int         precond_nbNonConstTWindows = -1;
  } _defaultValuesIQNILS;

  const struct DefaultValuesIMVJ {
    double      relaxationFactor           = 0.1;
    int         maxIterationsUsed          = 20;
    int         timeWindowsReused          = 0;
    int         filter                     = Acceleration::QR2FILTER;
    double      singularityLimit           = 1e-2;
    std::string preconditionerType         = "residual-sum";
    int         precond_nbNonConstTWindows = -1;
  } _defaultValuesIQNIMVJ;

  struct UserDefinitions {
    bool definedRelaxationFactor   = false;
    bool definedMaxIterationsUsed  = false;
    bool definedTimeWindowsReused  = false;
    bool definedFilter             = false;
    bool definedPreconditionerType = false;
  } _userDefinitions;

  void addTypeSpecificSubtags(xml::XMLTag &tag);
  void addCommonIQNSubtags(xml::XMLTag &tag);
};
} // namespace acceleration
} // namespace precice
