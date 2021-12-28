#pragma once

#include <Eigen/Core>

#include "DataContext.hpp"
#include "logging/Logger.hpp"
#include "time/SharedPointer.hpp"
#include "time/Time.hpp"

namespace precice {
namespace impl {

/**
 * @brief Stores one Data object with related mesh. Additionally stores Waveform associated with _providedData.
 *
 * Derived from DataContext
 */
class ReadDataContext : public DataContext {
public:
  ReadDataContext(
      mesh::PtrData data,
      mesh::PtrMesh mesh,
      int           interpolationOrder = time::Time::UNDEFINED_INTERPOLATION_ORDER);

  /**
   * @brief links a read mapping and the mesh context the read mapping requires to this data context
   *
   * @param[in] mappingContext provides context of read mapping
   * @param[in] fromMeshContext provides context of mesh this read mapping is mapping to (_toData)
   */
  void configureForReadMapping(MappingContext mappingContext, MeshContext fromMeshContext);

  /// helper function for mapping of data
  void mapReadData();

  /// helper for mapReadDataTo
  void mapReadDataTo();

  /**
   * @brief Allows to sample data at a given point in time insize of the time window
   *
   * @param normalizedDt defines point in time where waveform will be sampled. Must be normalized to [0,1], where 0 refers to the beginning and 1 to the end of the window.
   */
  Eigen::VectorXd sampleWaveformAt(double normalizedDt);

  /// for initializing waveforms of the context.
  void initializeWaveform();

  /// for moving to the next time window and updating waveform correspondingly.
  void moveToNextWindow();

private:
  logging::Logger _log{"impl::ReadDataContext"};

  /// Waveform wrapped by this ReadDataContext.
  time::PtrWaveform _providedWaveform;
};

} // namespace impl
} // namespace precice