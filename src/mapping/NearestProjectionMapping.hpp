#pragma once

#include "Mapping.hpp"
#include <list>
#include <vector>

namespace precice {
namespace mapping {

/**
 * @brief Mapping using orthogonal projection to nearest triangle/edge/vertex and
 *        linear interpolation from projected point.
 */
class NearestProjectionMapping : public Mapping
{
public:

  /**
   * @brief Constructor, taking mapping constraint.
   */
  NearestProjectionMapping ( Constraint constraint, int dimensions );

  /**
   * @brief Destructor, empty.
   */
  virtual ~NearestProjectionMapping() {}

  /**
   * @brief Computes the projections and interpolation relations.
   */
  virtual void computeMapping();

  virtual bool hasComputedMapping();

  /**
   * @brief Removes a computed mapping.
   */
  virtual void clear();

  /**
   * @brief Uses projection and interpolation relations to map data values.
   *
   * Preconditions:
   * - computeMapping() has been called
   *
   * @param inputDataID [IN] Data ID of input data values to be mapped from.
   * @param outputDataID [IN] Data ID of output data values to be mapped to.
   */
  virtual void map (
    int inputDataID,
    int outputDataID );

  virtual bool doesVertexContribute(int vertexID);
  virtual bool isProjectionMapping();

private:

  static tarch::logging::Log _log;

  typedef std::list<query::InterpolationElement> InterpolationElements;
  std::vector<InterpolationElements> _weights;

  bool _hasComputedMapping;
};

}} // namespace precice, mapping
