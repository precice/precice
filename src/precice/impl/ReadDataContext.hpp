#pragma once

#include <Eigen/Core>

#include "DataContext.hpp"
#include "logging/Logger.hpp"
#include "time/SharedPointer.hpp"
#include "time/Time.hpp"

namespace precice {
namespace impl {

/**
 * @brief Stores one Data object with related mesh. Context is used to be read from and potentially provides a read mapping. Additionally stores Waveform associated with _providedData.
 *
 * Derived from DataContext
 */
class ReadDataContext : public DataContext {
public:
  /**
   * @brief Construct a new ReadDataContext object without a mapping.
   *
   * @param data Data associated with this ReadDataContext.
   * @param mesh Mesh associated with this ReadDataContext.
   * @param interpolationOrder Order of the Waveform stored by this ReadDataContext.
   */
  ReadDataContext(
      mesh::PtrData data,
      mesh::PtrMesh mesh,
      int           interpolationOrder = time::Time::UNDEFINED_INTERPOLATION_ORDER);

  /**
   * @brief Links a MappingContext for a read mapping and the MeshContext the read mapping requires to this DataContext.
   *
   * A read mapping maps _fromData to _providedData. A ReadDataContext already has _providedData, but additionally requires _fromData.
   *
   * @param[in] mappingContext Context of read mapping
   * @param[in] fromMeshContext Context of mesh this read mapping is mapping from (_fromData)
   */
  void configureForReadMapping(MappingContext mappingContext, MeshContext fromMeshContext);

  /**
   * @brief Performs the mapping associated to this ReadDataContext. Called by SolverInterfaceImpl::mapReadData on all ReadDataContext objects.
   */
  void mapReadData();

  /**
   * @brief Performs the mapping associated to this ReadDataContext. Called by SolverInterfaceImpl::mapReadDataTo.
   */
  void mapReadDataTo();

  /**
   * @brief Allows to sample data at a given point in time insize of the time window
   *
   * @param normalizedDt defines point in time where waveform will be sampled. Must be normalized to [0,1], where 0 refers to the beginning and 1 to the end of the window.
   */
  Eigen::VectorXd sampleWaveformAt(double normalizedDt);

  /**
   * @brief Initializes the _providedWaveform as a constant function with values from _providedData.
   */
  void initializeWaveform();

  /**
   * @brief Updates _providedWaveform when moving to the next time window.
   */
  void moveToNextWindow();

private:
  logging::Logger _log{"impl::ReadDataContext"};

  /// Waveform wrapped by this ReadDataContext.
  time::PtrWaveform _providedWaveform;
};

} // namespace impl
} // namespace precice
