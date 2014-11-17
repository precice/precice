// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_MAPPING_NEARESTNEIGHBORMAPPING_HPP_
#define PRECICE_MAPPING_NEARESTNEIGHBORMAPPING_HPP_

#include "mapping/Mapping.hpp"
#include "tarch/logging/Log.h"
#include <vector>

namespace precice {
namespace mapping {

/**
 * @brief Mapping using nearest neighboring vertices.
 */
class NearestNeighborMapping : public Mapping
{
public:

  /**
   * @brief Constructor.
   *
   * @param constraint [IN] Specifies mapping to be consistent or conservative.
   */
  NearestNeighborMapping ( Constraint constraint );

  /**
   * @brief Destructor, empty.
   */
  virtual ~NearestNeighborMapping() {}

  /**
   * @brief Computes the mapping coefficients from the in- and output mesh.
   */
  virtual void computeMapping();

  /**
   * @brief Returns true, if computeMapping() has been called.
   */
  virtual bool hasComputedMapping();

  /**
   * @brief Removes a computed mapping.
   */
  virtual void clear();

  /**
   * @brief Maps input data to output data from input mesh to output mesh.
   */
  virtual void map (
    int inputDataID,
    int outputDataID );

  virtual bool doesVertexContribute(int vertexID);

private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  // @brief Flag to indicate whether computeMapping() has been called.
  bool _hasComputedMapping;

  // @brief Computed output vertex indices to map data from input vertices to.
  std::vector<int> _vertexIndices;
};

}} // namespace precice, mapping

#endif /* PRECICE_MAPPING_NEARESTNEIGHBORMAPPING_HPP_ */
