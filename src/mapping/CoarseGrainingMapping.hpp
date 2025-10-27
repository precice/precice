#pragma once

#include <string>
#include <vector>
#include "logging/Logger.hpp"
#include "mapping/Mapping.hpp"

namespace precice::mapping {

// Forward declaration
namespace impl {
class LucyKernelFunction;
}

/// Mapping using nearest neighboring vertices and (eventually) their local gradient values.
/// Base class for Nearest Neighbor Mapping and Nearest Neighbor Gradient
class CoarseGrainingMapping : public Mapping {
public:
  /**
   * @brief Constructor.
   *
   * @param[in] constraint Specifies mapping to be consistent or conservative.
   * @param[in] meshDimension Dimension of the meshes
   * @param[in] grainDimension Dimension of the coarse graining function
   */
  CoarseGrainingMapping(Constraint constraint, int meshDimension, int grainDimension, double functionRadius);

  /// Computes the mapping coefficients from the in- and output mesh.
  void computeMapping() final override;

  /// For reading data just-in-time (only consistent at the moment)
  void mapConsistentAt(const Eigen::Ref<const Eigen::MatrixXd> &coordinates, const impl::MappingDataCache &cache, Eigen::Ref<Eigen::MatrixXd> values) final override;

  /// For writing data just-in-time (only conservative at the moment)
  void mapConservativeAt(const Eigen::Ref<const Eigen::MatrixXd> &coordinates, const Eigen::Ref<const Eigen::MatrixXd> &source, impl::MappingDataCache &cache, Eigen::Ref<Eigen::MatrixXd> target) final override;

  /// Removes a computed mapping.
  void clear() final override;

  void tagMeshFirstRound() final override;

  void tagMeshSecondRound() final override;

  /// name of the coarse-graining mapping
  std::string getName() const final override;

private:
  mutable logging::Logger _log{"mapping::CoarseGraining"};

  /// @copydoc Mapping::mapConservative
  void mapConservative(const time::Sample &inData, Eigen::VectorXd &outData) override;

  /// @copydoc Mapping::mapConsistent
  void mapConsistent(const time::Sample &inData, Eigen::VectorXd &outData) override;

  std::unique_ptr<impl::LucyKernelFunction> _lucyFunction;
};
} // namespace precice::mapping
