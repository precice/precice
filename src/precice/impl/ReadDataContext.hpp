#pragma once

#include <Eigen/Core>

#include "DataContext.hpp"
#include "cplscheme/ImplicitData.hpp"
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
   * @brief Gets degree of waveform
   *
   * @return int degree of waveform
   */
  int getWaveformDegree() const;

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
   * @param[in] time Point in time where waveform is sampled.
   * @param[in] values read data associated with given indices for time \ref time will be returned into this span
   */
  void readValues(::precice::span<const VertexID> vertices, double time, ::precice::span<double> values) const;

  /// Disable copy construction
  ReadDataContext(const ReadDataContext &copy) = delete;

  /// Disable assignment construction
  ReadDataContext &operator=(const ReadDataContext &assign) = delete;

  /// Move constructor, use the implicitly declared.
  ReadDataContext(ReadDataContext &&) = default;
  ReadDataContext &operator=(ReadDataContext &&) = default;

  /**
   * @brief Removes all toData samples from mappings
   */
  void clearToDataFor(const cplscheme::ImplicitData &from);

  /**
   * @brief Trims all toData of associated mappings after the given t
   *
   * @param[in] t the time after which to trim data
   */
  void trimToDataAfterFor(const cplscheme::ImplicitData &from, double t);

private:
  static logging::Logger _log;
};

} // namespace impl
} // namespace precice
