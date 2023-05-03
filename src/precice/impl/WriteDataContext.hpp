#pragma once

#include <Eigen/Core>

#include "DataContext.hpp"
#include "logging/Logger.hpp"
#include "time/Sample.hpp"

namespace precice {
namespace impl {

/**
 * @brief Stores one Data object with related mesh. Context stores data to be written to and potentially provides a write mapping.
 *
 * Derived from DataContext
 */
class WriteDataContext : public DataContext {
public:
  /**
   * @brief Construct a new WriteDataContext object without a mapping.
   *
   * @param data Data associated with this WriteDataContext.
   * @param mesh Mesh associated with this WriteDataContext.
   */
  WriteDataContext(
      mesh::PtrData data,
      mesh::PtrMesh mesh);

  /**
   * @brief Resets provided data, writeDataBuffer, (if mapping exists) fromData or toData, and (optionally) storage.
   *
   * @param atEndOfWindow if true, also the Storage will be reset (useful at end of window to trim storage).
   * @param isTimeWindowComplete if true, overwrite sample at front of Storage with sample at back (basically a Storage::move with constant).
   */
  void resetData(bool atEndOfWindow, bool isTimeWindowComplete);

  /**
   * @brief Store values in _writeDataBuffer
   *
   * @param[in] vertices ids of data
   * @param[in] values values of data
   */
  void writeValuesIntoDataBuffer(::precice::span<const VertexID> vertices, ::precice::span<const double> values);

  /**
   * @brief Store gradients in _writeDataBuffer
   *
   * @param[in] vertices ids of data
   * @param[in] gradients gradients of data
   */
  void writeGradientsIntoDataBuffer(::precice::span<const VertexID> vertices, ::precice::span<const double> gradients);

  void resizeBufferTo(int size);

  /**
   * @brief Store data from _writeDataBuffer in persistent storage
   *
   * @param[in] currentTime time data should be associated with
   */
  void storeBufferedData(double currentTime);

  /**
   * @brief Adds a MappingContext and the MeshContext required by the write mapping to the corresponding WriteDataContext data structures.
   *
   * A write mapping maps _providedData to _toData. A WriteDataContext already has _providedData, but additionally requires _toData.
   *
   * @param[in] mappingContext Context of write mapping
   * @param[in] meshContext Context of mesh this write mapping is mapping to (_toData)
   */
  void appendMappingConfiguration(MappingContext &mappingContext, const MeshContext &meshContext);

  /// Disable copy construction
  WriteDataContext(const WriteDataContext &copy) = delete;

  /// Disable assignment construction
  WriteDataContext &operator=(const WriteDataContext &assign) = delete;

  /// Move constructor, use the implicitly declared.
  WriteDataContext(WriteDataContext &&) = default;
  WriteDataContext &operator=(WriteDataContext &&) = default;

private:
  static logging::Logger _log;

  /// @brief Buffer to store written data until it is copied to _providedData->timeStepsStorage()
  time::Sample _writeDataBuffer;
};

} // namespace impl
} // namespace precice
