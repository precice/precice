#pragma once

#include "action/Action.hpp"
#include "action/SharedPointer.hpp"
#include "xml/XMLTag.hpp"
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"
#include <string>
#include <list>

namespace precice {
namespace action {

/**
 * @brief Configures an Action subclass object.
 */
class ActionConfiguration : public xml::XMLTag::Listener
{
public:

  ActionConfiguration (
    xml::XMLTag&                    parent,
    const mesh::PtrMeshConfiguration& meshConfig );

  /**
   * @brief Callback function required for use of automatic configuration.
   *
   * @return True, if successful.
   */
  virtual void xmlTagCallback ( xml::XMLTag& callingTag );

  /**
   * @brief Callback function required for use of automatic configuration.
   *
   * @return True, if successful.
   */
  virtual void xmlEndTagCallback ( xml::XMLTag& callingTag );

  /**
   * @brief Returns the id of the mesh used in the data action.
   */
  int getUsedMeshID() const;

  /**
   * @brief Returns the configured action.
   */
  const std::list<PtrAction>& actions() const
  {
    return _actions;
  }

  void resetActions()
  {
    _actions.clear();
  }

private:

  /**
   * @brief Stores configuration information temporarily to create the Action.
   */
  struct ConfiguredAction
  {
    std::string type;
    std::string timing;
    std::string sourceData;
    std::string targetData;
    std::string mesh;
    double convergenceTolerance = 0;
    int maxIterations = 0;
    std::string path;
    std::string module;
  };

  mutable logging::Logger _log{"config::ActionConfiguration"};

  const std::string TAG = "action";

  const std::string NAME_DIVIDE_BY_AREA;
  const std::string NAME_MULTIPLY_BY_AREA;
  const std::string NAME_SCALE_BY_COMPUTED_DT_RATIO;
  const std::string NAME_SCALE_BY_COMPUTED_DT_PART_RATIO;
  const std::string NAME_SCALE_BY_DT;
  const std::string NAME_ADD_TO_COORDINATES;
  const std::string NAME_SUBTRACT_FROM_COORDINATES;
  const std::string NAME_COMPUTE_CURVATURE;
  const std::string NAME_PYTHON;

  const std::string TAG_SOURCE_DATA;
  const std::string TAG_TARGET_DATA;
  const std::string TAG_CONVERGENCE_TOLERANCE;
  const std::string TAG_MAX_ITERATIONS;
  const std::string TAG_MODULE_PATH;
  const std::string TAG_MODULE_NAME;

  const std::string ATTR_TYPE = "type";
  const std::string ATTR_TIMING = "timing";
  const std::string ATTR_NAME = "name";
  const std::string ATTR_VALUE = "value";
  const std::string ATTR_MESH = "mesh";

  const std::string VALUE_REGULAR_PRIOR;
  const std::string VALUE_REGULAR_POST;
  const std::string VALUE_ON_EXCHANGE_PRIOR;
  const std::string VALUE_ON_EXCHANGE_POST;
  const std::string VALUE_ON_TIMESTEP_COMPLETE_POST;

  mesh::PtrMeshConfiguration _meshConfig;

  ConfiguredAction _configuredAction;

  std::list<PtrAction> _actions;

//  /**
//   * @brief Adds all required subtags to the main action tag.
//   */
//  void addSubtags (
//    std::list<xml::XMLTag>& tags,
//    const std::string&        type );

  /**
   * @brief Creates the Action object.
   */
  void createAction();

  Action::Timing getTiming() const;
};

}} // namespace precice, action
