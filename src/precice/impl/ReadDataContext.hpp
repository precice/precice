#pragma once

#include <Eigen/Core>

#include "DataContext.hpp"
#include "logging/Logger.hpp"

namespace precice {
namespace impl {

/**
 * @brief Stores one Data object with related mesh. Context stores data to be read from and potentially provides a read mapping. Additionally stores Waveform object associated with _providedData.
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
   */
  ReadDataContext(
      mesh::PtrData data,
      mesh::PtrMesh mesh);

  /**
   * @brief Gets _interpolationOrder of _waveform
   *
   * @return _interpolationOrder of _waveform
   */
  int getInterpolationOrder() const;

  /**
   * @brief Adds a MappingContext and the MeshContext required by the read mapping to the corresponding ReadDataContext data structures.
   *
   * A read mapping maps _fromData to _providedData. A ReadDataContext already has _providedData, but additionally requires _fromData.
   *
   * @param[in] mappingContext Context of read mapping
   * @param[in] meshContext Context of mesh this read mapping is mapping from (_fromData)
   *
   * @note Note that read mapping contexts have to be unique and and this function may only be called once for a given ReadDataContext
   */
  void appendMappingConfiguration(MappingContext &mappingContext, const MeshContext &meshContext) override;

  /**
   * @brief Samples data at a given point in time within the current time window for given indices
   *
   * @param[in] vertices vertex ids
   * @param[in] normalizedDt Point in time where waveform is sampled. Must be normalized to [0,1], where 0 refers to the beginning and 1 to the end of the current time window.
   * @param[in] values read data associated with given indices for time normalizedDt will be returned into this span
   */
  void readValues(::precice::span<const VertexID> vertices, double normalizedDt, ::precice::span<double> values) const;

  /// Disable copy construction
  ReadDataContext(const ReadDataContext &copy) = delete;

  /// Disable assignment construction
  ReadDataContext &operator=(const ReadDataContext &assign) = delete;

  /// Move constructor, use the implicitly declared.
  ReadDataContext(ReadDataContext &&) = default;
  ReadDataContext &operator=(ReadDataContext &&) = default;

private:
  static logging::Logger _log;
};

} // namespace impl
} // namespace precice
