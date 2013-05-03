// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CONFIG_ACTIONCONFIGURATION_HPP_
#define PRECICE_CONFIG_ACTIONCONFIGURATION_HPP_

#include "action/Action.hpp"
#include "action/SharedPointer.hpp"
#include "utils/xml/XMLTag.hpp"
#include "tarch/logging/Log.h"
#include "mesh/SharedPointer.hpp"
#include <string>
#include <list>

namespace precice {
namespace action {

/**
 * @brief Configures an Action subclass object.
 */
class ActionConfiguration : public utils::XMLTag::Listener
{
public:

  /**
   * @brief Returns the name of the enclosing XML-tag.
   */
  //static const std::string& getTag();

  /**
   * @brief Constructor.
   */
  ActionConfiguration (
    utils::XMLTag&                    parent,
    const mesh::PtrMeshConfiguration& meshConfig );

//  /**
//   * @brief Reads the information parsed from an xml-file.
//   */
//  bool parseSubtag ( utils::XMLTag::XMLReader* xmlReader );

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
   * @brief Returns true, if configuration has validly taken place.
   */
  //bool isValid() const;

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
    double convergenceTolerance;
    int maxIterations;
    std::string path;
    std::string module;

    ConfiguredAction ()
    : type (), timing(), sourceData(), targetData(), mesh(),
      convergenceTolerance(0.0), maxIterations(0), path(), module()
    {}
  };

  static tarch::logging::Log _log;

  const std::string TAG;

  const std::string NAME_DIVIDE_BY_AREA;
  const std::string NAME_MULTIPLY_BY_AREA;
  const std::string NAME_SCALE_BY_COMPUTED_DT_RATIO;
  const std::string NAME_SCALE_BY_COMPUTED_DT_PART_RATIO;
  const std::string NAME_SCALE_BY_DT;
  const std::string NAME_ADD_TO_COORDINATES;
  const std::string NAME_SUBTRACT_FROM_COORDINATES;
  const std::string NAME_COMPUTE_CURVATURE;
  const std::string NAME_BALANCE_VERTEX_POSITIONS;
  const std::string NAME_PYTHON;

  const std::string TAG_SOURCE_DATA;
  const std::string TAG_TARGET_DATA;
  const std::string TAG_CONVERGENCE_TOLERANCE;
  const std::string TAG_MAX_ITERATIONS;
  const std::string TAG_MODULE_PATH;
  const std::string TAG_MODULE_NAME;

  const std::string ATTR_TYPE;
  const std::string ATTR_TIMING;
  const std::string ATTR_NAME;
  const std::string ATTR_VALUE;
  const std::string ATTR_MESH;

  const std::string VALUE_REGULAR_PRIOR;
  const std::string VALUE_REGULAR_POST;
  const std::string VALUE_ON_EXCHANGE_PRIOR;
  const std::string VALUE_ON_EXCHANGE_POST;
  const std::string VALUE_ON_TIMESTEP_COMPLETE_POST;

  //bool _isValid;

  mesh::PtrMeshConfiguration _meshConfig;

  ConfiguredAction _configuredAction;

  std::list<PtrAction> _actions;

//  /**
//   * @brief Adds all required subtags to the main action tag.
//   */
//  void addSubtags (
//    std::list<utils::XMLTag>& tags,
//    const std::string&        type );

  /**
   * @brief Creates the Action object.
   */
  void createAction();

  Action::Timing getTiming() const;
};

}} // namespace precice, action

#endif /* PRECICE_CONFIG_ACTIONCONFIGURATION_HPP_ */
