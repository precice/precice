// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CONFIG_ACTIONCONFIGURATION_HPP_
#define PRECICE_CONFIG_ACTIONCONFIGURATION_HPP_

#include "precice/impl/SharedPointer.hpp"
#include "utils/xml/XMLTag.hpp"
#include "tarch/irr/XML.h"
#include "tarch/logging/Log.h"
#include "utils/PointerVector.hpp"
#include "precice/impl/MeshContext.hpp"
#include "precice/impl/AbstractDataAction.hpp"
#include <string>

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace config {

/**
 * @brief Configures an AbstractDataAction subclass object.
 */
class ActionConfiguration
{
public:

  /**
   * @brief Returns the name of the enclosing XML-tag.
   */
  static const std::string & getTag ();

  /**
   * @brief Constructor.
   */
  ActionConfiguration ( utils::ptr_vector<impl::MeshContext> & usedMeshes );

  /**
   * @brief Reads the information parsed from an xml-file.
   */
  bool parseSubtag ( tarch::irr::io::IrrXMLReader * xmlReader );

  /**
   * @brief Callback function required for use of automatic configuration.
   *
   * @return True, if successful.
   */
  void xmlTagCallback (
    utils::XMLTag         & callingTag,
    tarch::irr::io::IrrXMLReader * xmlReader );

  /**
   * @brief Callback function required for use of automatic configuration.
   *
   * @return True, if successful.
   */
  void xmlEndTagCallback (
    utils::XMLTag &                callingTag,
    tarch::irr::io::IrrXMLReader * xmlReader );

  /**
   * @brief Returns true, if configuration has validly taken place.
   */
  bool isValid () const;

  /**
   * @brief Returns the configured action.
   */
  const impl::PtrAbstractDataAction & getAction () const;

private:

  static tarch::logging::Log _log;

  const std::string TAG_SOURCE_DATA;
  const std::string TAG_TARGET_DATA;
  const std::string TAG_MESH;
  const std::string TAG_CONVERGENCE_TOLERANCE;
  const std::string TAG_MAX_ITERATIONS;

  const std::string ATTR_TYPE;
  const std::string ATTR_TIMING;
  const std::string ATTR_NAME;
  const std::string ATTR_VALUE;

  const std::string VALUE_REGULAR_PRIOR;
  const std::string VALUE_REGULAR_POST;
  const std::string VALUE_ON_EXCHANGE_PRIOR;
  const std::string VALUE_ON_EXCHANGE_POST;

  const std::string VALUE_DIVIDE_BY_AREA;
  const std::string VALUE_MULTIPLY_BY_AREA;
  const std::string VALUE_SCALE_BY_COMPUTED_DT;
  const std::string VALUE_SCALE_BY_COMPUTED_DT_PART;
  //const std::string VALUE_SET_AS_COORDINATES;
  const std::string VALUE_ADD_TO_COORDINATES;
  const std::string VALUE_SUBTRACT_FROM_COORDINATES;
  const std::string VALUE_COMPUTE_CURVATURE;
  const std::string VALUE_BALANCE_VERTEX_POSITIONS;

  bool _isValid;

  utils::ptr_vector<impl::MeshContext> & _usedMeshes;

  impl::PtrAbstractDataAction _action;

  /**
   * @brief Stores configuration information temporarily to create the Action.
   */
  struct Configured {
    std::string type;
    std::string timing;
    std::string sourceData;
    std::string targetData;
    std::string mesh;
    double convergenceTolerance;
    int maxIterations;

    Configured ()
    : type (), timing(), sourceData(), targetData(), mesh(),
      convergenceTolerance(0.0), maxIterations(0)
    {}
  } _configured;

  /**
   * @brief Adds all required subtags to the main action tag.
   */
  void addSubtags (
    utils::XMLTag &     callingTag,
    const std::string & type );

  /**
   * @brief Creates the Action object.
   */
  void createAction ();

  impl::AbstractDataAction::TimingConstants getTiming () const;
};

}} // namespace precice, config

#endif /* PRECICE_CONFIG_ACTIONCONFIGURATION_HPP_ */
