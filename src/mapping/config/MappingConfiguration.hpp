// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_MAPPING_MAPPINGCONFIGURATION_HPP_
#define PRECICE_MAPPING_MAPPINGCONFIGURATION_HPP_

#include "mapping/SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"
#include "tarch/logging/Log.h"
#include "utils/xml/XMLTag.hpp"
#include <string>
#include <vector>

namespace precice {
namespace mapping {

/**
 * @brief Performs XML configuration and holds configured mappings.
 */
class MappingConfiguration : public utils::XMLTag::Listener
{
public:

  /**
   * @brief Constants defining the direction of a mapping.
   */
  enum Direction
  {
    WRITE,
    READ
  };

  enum Timing
  {
    INITIAL,
    ON_ADVANCE,
    ON_DEMAND,
    INCREMENTAL
  };

  /**
   * @brief Configuration data for one mapping.
   */
  struct ConfiguredMapping
  {
    // @brief Mapping object.
    PtrMapping mapping;
    // @brief Remote mesh to map to/from
    mesh::PtrMesh mesh;
    // @brief Direction of mapping (important to set input and output mesh).
    Direction direction;
    // @brief When the mapping should be executed.
    Timing timing;
    // @brief Is true, if the mapping is done incremental for each value.
    //bool isIncremental;
  };

  // @brief Name of xml tag for this class in configuration file
  //static const std::string& getTag();

  /**
   * @brief Constructor.
   */
  MappingConfiguration (
    utils::XMLTag&                    parent,
    const mesh::PtrMeshConfiguration& meshConfiguration );

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
   * @returns Returns true, if the xml-file parsing was successful.
   */
  //bool isValid() const;

  /**
   * @brief Returns all configured mappings.
   */
  const std::vector<ConfiguredMapping>& mappings();

  /**
   * @brief Adds a mapping to the configuration.
   */
  void addMapping (
    const PtrMapping&    mapping,
    const mesh::PtrMesh& mesh,
    Direction            direction,
    //bool                 isIncremental,
    Timing               timing );

  void resetMappings() { _mappings.clear(); }

private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  const std::string TAG;

  const std::string ATTR_DIRECTION;
  const std::string ATTR_MESH;
  const std::string ATTR_TIMING;
  //const std::string ATTR_INCREMENTAL;
  const std::string ATTR_TYPE;
  const std::string ATTR_CONSTRAINT;
  const std::string ATTR_SHAPE_PARAM;
  const std::string ATTR_SUPPORT_RADIUS;

  const std::string VALUE_WRITE;
  const std::string VALUE_READ;
  const std::string VALUE_CONSISTENT;
  const std::string VALUE_CONSERVATIVE;
  const std::string VALUE_NEAREST_NEIGHBOR;
  const std::string VALUE_NEAREST_PROJECTION;
  const std::string VALUE_RBF_TPS;
  const std::string VALUE_RBF_MULTIQUADRICS;
  const std::string VALUE_RBF_INV_MULTIQUADRICS;
  const std::string VALUE_RBF_VOLUME_SPLINES;
  const std::string VALUE_RBF_GAUSSIAN;
  const std::string VALUE_RBF_CTPS_C2;
  const std::string VALUE_RBF_CPOLYNOMIAL_C0;
  const std::string VALUE_RBF_CPOLYNOMIAL_C6;
  const std::string VALUE_TIMING_INITIAL;
  const std::string VALUE_TIMING_ON_ADVANCE;
  const std::string VALUE_TIMING_ON_DEMAND;
  const std::string VALUE_TIMING_INCREMENTAL;

  mesh::PtrMeshConfiguration _meshConfig;

  //bool _isValid;

  std::vector<ConfiguredMapping> _mappings;

  ConfiguredMapping createMapping (
    const std::string& direction,
    const std::string& type,
    const std::string& constraint,
    const std::string& meshName,
    Timing             timing,
    //bool               incremental,
    double             shapeParameter,
    double             supportRadius ) const;

  void checkDuplicates ( const ConfiguredMapping& mapping );

  Timing getTiming(const std::string& timing) const;
};

}} // namespace mapping, config

#endif /* PRECICE_MAPPING_MAPPINGCONFIGURATION_HPP_ */
