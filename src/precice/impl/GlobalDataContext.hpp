#pragma once

#include <string>
// #include "MappingContext.hpp"
// #include "MeshContext.hpp"
#include <Eigen/Core>
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"
#include "time/SharedPointer.hpp"
#include "time/Time.hpp"

namespace precice {
namespace impl {

/**
 * @brief Stores one Global (meshless) Data object.
 *
 * - GlobalDataContext is basically a stripped-down equivalent of DataContext.
 *   Unlike DataContext, it has no associated mappings or meshes.
 */
class GlobalDataContext {

public:
  /**
   * @brief Construct a new GlobalDataContext with time interpolation.
   *
   * @param data Data associated with this DataContext.
   */
  GlobalDataContext(mesh::PtrGlobalData data);

  /**
   * @brief Get the Name of _providedData.
   *
   * @return std::string Name of _providedData.
   */
  std::string getDataName() const;

  /**
   * @brief Get the dimensions of _providedData.
   *
   * @return int Dimensions of _providedData.
   */
  int getDataDimensions() const;

  /**
   * @brief Resets provided data and (if mapping exists) fromData or toData.
   */
  void resetData();

  /**
   * @brief Get _providedData member.
   *
   * @return mesh::PtrGlobalData _providedData.
   */
  mesh::PtrGlobalData providedData() const;

  /**
   * @brief Gets _interpolationOrder of _waveform
   *
   * @return _interpolationOrder of _waveform
   */
  int getInterpolationOrder() const;

  /**
   * @brief Samples data at a given point in time within the current time window
   *
   * @param normalizedDt Point in time where waveform is sampled. Must be normalized to [0,1], where 0 refers to the beginning and 1 to the end of the current time window.
   */
  Eigen::VectorXd sampleWaveformAt(double normalizedDt);

  /**
   * @brief Initializes the _waveform as a constant function with values from _providedData.
   */
  void initializeWaveform();

  /**
   * @brief Updates _waveform when moving to the next time window.
   */
  void moveToNextWindow();

  /**
   * @brief Return the direction (read or write) for this GlobalDataContext.
   */
  std::string getDirection();

  /**
   * @brief Stores _providedData as first sample of _waveform.
   */
  void storeDataInWaveform();

private:
  static logging::Logger _log;

  /// Unique data this context is associated with
  mesh::PtrGlobalData _providedData;

  /// Whether this is "read" context or "write" context
  std::string _direction;

  /// Waveform wrapped by this GlobalDataContext.
  time::PtrWaveform _waveform;
};

} // namespace impl
} // namespace precice
