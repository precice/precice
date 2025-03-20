#pragma once

#include <Eigen/Core>

#include "DataContext.hpp"
#include "logging/Logger.hpp"
#include "time/Sample.hpp"

namespace precice::impl {

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

  /// Resets the writeDataBuffer and the mapping data cache
  void resetBufferedData();

  /**
   * @brief Removes stample before \ref time and (if mapping exists) fromData or toData
   *
   * @param time the point in time after which to remove samples
   */
  void trimAfter(double time);

  /**
   * @brief Forwards the just-in-time mapping API call for writing data to the data context
   *
   * writeAndMapValues then forwards the arguments along with a target buffer
   * to the \p justInTimeMapping this DataContext holds. The function always forwards both,
   * the \p _writeDataBuffer and the \p mappingCache. The \p justInTimeMapping then decides,
   * whether the data can directly be stored in the actual target \p _writeDataBuffer or whether
   * intermediate results are stored in the \p mappingCache.
   *
   * If written data is stored in the \p mappingCache, the function \p completeJustInTimeMapping
   * is responsible for evaluating the cache and storing the result in the \p _writeDataBuffer.
   *
   * See also the documentation of the impl::MappingDataCache class.
   *
   * @param[in] coordinates Provided by the user through the API function
   * @param[in] values Containing the write data values given by the user
   */
  void writeAndMapValues(::precice::span<const double> coordinates, ::precice::span<const double> values);

  /**
   * @brief Evaluates the MappingDataCache and stores the result in the \p _writeDataBuffer
   *
   * If a just-in-time mapping stores intermediate results in the \p mappingCache, we need
   * to finalize the mapping in the cache and store the mapped values in the \p _writeDataBuffer.
   *
   * See also \p writeAndMapValues and the documentation of Mapping::completeJustInTimeMapping as well as
   * the impl::MappingDataCache class for more details.
   *
   * This function needs to be called before calling \p storeBufferedData, which stores the buffer
   * into the time step storage.
   */
  void completeJustInTimeMapping();

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

  void resizeBufferTo(int size, bool invalidateBufferedData = false);

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
  void appendMappingConfiguration(MappingContext &mappingContext, const MeshContext &meshContext) override;

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

} // namespace precice::impl
