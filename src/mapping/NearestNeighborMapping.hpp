#pragma once

#include <string>
#include <vector>
#include "logging/Logger.hpp"
#include "mapping/NearestNeighborBaseMapping.hpp"

namespace precice::mapping {

/// Mapping using nearest neighboring vertices.
class NearestNeighborMapping : public NearestNeighborBaseMapping {
public:
  /**
   * @brief Constructor.
   *
   * @param[in] constraint Specifies mapping to be consistent or conservative.
   * @param[in] dimensions Dimensionality of the meshes
   */
  NearestNeighborMapping(Constraint constraint, int dimensions);

  /// name of the nn mapping
  std::string getName() const final override;

  /// For reading data just-in-time (only consistent at the moment)
  void mapConsistentAt(const Eigen::Ref<const Eigen::MatrixXd> &coordinates, const MappingDataCache &cache, Eigen::Ref<Eigen::MatrixXd> values) final override;

  /// For writing data just-in-time (only conservative at the moment)
  void mapConservativeAt(const Eigen::Ref<const Eigen::MatrixXd> &coordinates, const Eigen::Ref<const Eigen::MatrixXd> &source, Eigen::Ref<Eigen::MatrixXd> target);

protected:
  /// @copydoc Mapping::mapConservative
  void mapConservative(const time::Sample &inData, Eigen::VectorXd &outData) final override;

  /// @copydoc Mapping::mapConsistent
  void mapConsistent(const time::Sample &inData, Eigen::VectorXd &outData) final override;
};

} // namespace precice::mapping
