#pragma once

#include <string>
#include "MappingContext.hpp"
#include "MeshContext.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {

namespace testing {
// Forward declaration to friend the boost test struct
class DataContextFixture;
} // namespace testing

namespace impl {

/**
 * @brief Stores one Data object with related mesh.
 *
 * - If this dataContext is associated with a mapping, fromData and toData will be set correspondingly.
 *   One of the two must be equal to providedData. fromData and toData must be different.
 * - If this dataContext is not associated with a mapping, fromData and toData will be unset.
 */
class DataContext {
  friend class testing::DataContextFixture; // Make the fixture friend of this class
public:
  DataContext(mesh::PtrData data, mesh::PtrMesh mesh);

  mesh::PtrData providedData();

  std::string getDataName() const;

  int getProvidedDataID() const;

  int getFromDataID() const;

  /// for copying read data into waveform, if no mapping exists
  void moveProvidedDataToProvidedWaveform();

  int getToDataID() const;

  std::string getMeshName() const;

  int getMeshID() const;

  /**
   * @brief links a read mapping and the mesh context the read mapping requires to this data context
   *
   * @param[in] mappingContext provides context of read mapping
   * @param[in] fromMeshContext provides context of mesh this read mapping is mapping to (_toData)
   */
  void configureForReadMapping(MappingContext mappingContext, MeshContext fromMeshContext);

  /**
   * @brief links a write mapping and the mesh context the write mapping requires to this data context
   *
   * @param[in] mappingContext provides context of write mapping
   * @param[in] fromMeshContext provides context of mesh this write mapping is mapping from (_fromData)
   */
  void configureForWriteMapping(MappingContext mappingContext, MeshContext toMeshContext);

  bool hasMapping() const;

  const MappingContext mappingContext() const;

  /// for updating waveforms after communication
  void initializeContextWaveforms();

  /// for communication of read and write data
  void sampleWaveformInToData();

  void storeFromDataInWaveform();

  /// shift data in time
  void moveToNextWindow();

  /// helper function to map a Waveform sample
  void mapWaveformSample();

  /**
   * @brief Allows to sample data at a given point in time insize of the time window
   * 
   * @param normalizedDt defines point in time where waveform will be sampled. Must be normalized to [0,1], where 0 refers to the beginning and 1 to the end of the window.
   */
  Eigen::VectorXd sampleAt(double normalizedDt);

private:
  mutable logging::Logger _log{"impl::DataContext"};

  mesh::PtrMesh _mesh;

  // data this participant will write to and read from
  mesh::PtrData _providedData;

  mesh::PtrData _fromData;

  mesh::PtrData _toData;

  MappingContext _mappingContext;

  void resetProvidedData();

  void resetToData();

  bool hasReadMapping() const;

  bool hasWriteMapping() const;

  /// helper function for initializing waveforms from given data
  void initializeWaveform(mesh::PtrData initializingData);

  /// helper function for creating read and write mappings
  void setMapping(MappingContext mappingContext, mesh::PtrData fromData, mesh::PtrData toData);
};

} // namespace impl
} // namespace precice
