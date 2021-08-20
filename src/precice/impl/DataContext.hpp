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

  std::string getDataName() const;

  int getProvidedDataID() const;

  int getFromDataID() const;

  void resetProvidedData();

  void resetToData();

  int getToDataID() const;

  std::string getMeshName() const;

  int getMeshID() const;

  // links a read mapping to this data context
  void configureForReadMapping(MappingContext mappingContext, MeshContext meshContext);

  // links a write mapping to this data context
  void configureForWriteMapping(MappingContext mappingContext, MeshContext meshContext);

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

  int numberOfSamplesInWaveform();

  Eigen::VectorXd sampleAt(double dt, int timeWindows);

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
