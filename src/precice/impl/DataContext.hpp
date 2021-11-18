#pragma once

#include <string>
#include "MappingContext.hpp"
#include "MeshContext.hpp"
#include "mesh/SharedPointer.hpp"
#include "time/SharedPointer.hpp"

namespace precice {
namespace impl {

/**
 * @brief Stores one Data object with related mesh.
 *
 * - If this dataContext is associated with a mapping, fromData and toData will be set correspondingly.
 *   One of the two must be equal to providedData. fromData and toData must be different.
 * - If this dataContext is not associated with a mapping, fromData and toData will be unset.
 */
class DataContext {

public:
  DataContext(mesh::PtrData data, mesh::PtrMesh mesh);

  mesh::PtrData providedData();

  mesh::PtrData toData();

  std::string getDataName() const;

  int getProvidedDataID() const;

  int getFromDataID() const;

  void resetProvidedData();

  void resetToData();

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

  bool hasReadMapping() const;

  bool hasWriteMapping() const;

  const MappingContext mappingContext() const;

  // for updating waveforms after communication
  void initializeProvidedWaveform();

  void initializeFromWaveform();

  void initializeToWaveform();

  // for mapping waveforms
  void moveWaveformSampleToData(int sampleID);

  void moveDataToWaveformSample(int sampleID);

  // for communication of read and write data
  void sampleWaveformInToData();

  void storeFromDataInWaveform();

  // for actions
  void sampleWaveformInProvidedData();

  void storeProvidedDataInWaveform();

  // for copying read data into waveform
  void moveProvidedDataToProvidedWaveformSample(int sampleID);

  // shift data in time
  void moveProvidedWaveform();

  int sizeOfSampleStorageInWaveform();

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
  time::PtrWaveform _providedWaveform;

  mesh::PtrData _providedData;

  time::PtrWaveform _fromWaveform;

  mesh::PtrData _fromData;

  time::PtrWaveform _toWaveform;

  mesh::PtrData _toData;

  MappingContext _mappingContext;

  // helper function for initializing waveforms from given data
  void initializeWaveform(mesh::PtrData initializingData, time::PtrWaveform initializedWaveform);

  // helper functions for communication
  void sampleWaveformIntoData(mesh::PtrData targetData, time::PtrWaveform sourceWaveform, int sampleID = 0);

  void storeDataInWaveform(mesh::PtrData sourceData, time::PtrWaveform targetWaveform, int sampleID = 0);

  // helper function for creating read and write mappings
  void setMapping(MappingContext mappingContext, mesh::PtrData fromData, mesh::PtrData toData, time::PtrWaveform fromWaveform, time::PtrWaveform toWaveform);
};

} // namespace impl
} // namespace precice
