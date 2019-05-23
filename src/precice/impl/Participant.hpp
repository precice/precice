#pragma once

#include "precice/Constants.hpp"
#include "SharedPointer.hpp"
#include "action/SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"
#include "mapping/SharedPointer.hpp"
#include "io/config/ExportConfiguration.hpp"
#include "io/ExportContext.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "utils/PointerVector.hpp"
#include "partition/ReceivedPartition.hpp"
#include "utils/MasterSlave.hpp"
#include <string>

namespace precice {
  namespace impl {
    struct DataContext;
    struct MeshContext;
    struct MappingContext;
  }
}

// Forward declaration to friend the boost test struct
namespace PreciceTests {
  namespace Serial {
    struct TestConfiguration;
  }
}


namespace precice {
namespace impl {

/// Holds coupling state of one participating solver in coupled simulation.
class Participant
{
public:

  enum MappingConstants {
    MAPPING_LINEAR_CONSERVATIVE,
    MAPPING_LINEAR_CONSISTENT,
    MAPPING_DIRECT
  };

  static void resetParticipantCount();

  /**
   * @brief Constructor.
   *
   * @param[in] name Name of the participant. Has to be unique.
   */
  Participant (
    std::string                 name,
    mesh::PtrMeshConfiguration& meshConfig );

  virtual ~Participant();

  /// Returns the name of the participant.
  const std::string& getName() const;

  int getID() const;

  void addWriteData (
    const mesh::PtrData& data,
    const mesh::PtrMesh& mesh );

  void addReadData (
    const mesh::PtrData& data,
    const mesh::PtrMesh& mesh );

  const DataContext& dataContext ( int dataID ) const;

  DataContext& dataContext ( int dataID );

  const utils::ptr_vector<DataContext>& writeDataContexts() const;

  utils::ptr_vector<DataContext>& writeDataContexts();

  const utils::ptr_vector<DataContext>& readDataContexts() const;

  utils::ptr_vector<DataContext>& readDataContexts();

  bool isMeshUsed ( int meshID ) const;

  bool isDataUsed ( int dataID ) const;

  bool isDataRead ( int dataID ) const;

  bool isDataWrite ( int dataID ) const;

  const MeshContext& meshContext ( int meshID ) const;

  MeshContext& meshContext ( int meshID );

  const std::vector<MeshContext*>& usedMeshContexts() const;

  std::vector<MeshContext*>& usedMeshContexts();

  void addReadMappingContext(MappingContext* mappingContext);

  void addWriteMappingContext(MappingContext* mappingContext);

  const utils::ptr_vector<MappingContext>& readMappingContexts() const;

  const utils::ptr_vector<MappingContext>& writeMappingContexts() const;

  void addWatchPoint ( const PtrWatchPoint& watchPoint );

  std::vector<PtrWatchPoint>& watchPoints();

  /// Adds a mesh to be used by the participant.
  void useMesh (
    const mesh::PtrMesh&                          mesh,
    const Eigen::VectorXd&                        localOffset,
    bool                                          remote,
    const std::string&                            fromParticipant,
    double                                        safetyFactor,
    bool                                          provideMesh,
    partition::ReceivedPartition::GeometricFilter geoFilter);

  void addAction ( const action::PtrAction& action );

  std::vector<action::PtrAction>& actions();

  const std::vector<action::PtrAction>& actions() const;

  /// Adds an export context to export meshes and data.
  void addExportContext ( const io::ExportContext& context );

  /// Returns all export contexts for exporting meshes and data.
  const std::vector<io::ExportContext>& exportContexts() const;

  /// Returns true, if the participant uses a precice in form of a server.
  bool useServer();

  /// Sets the client-server com. for the participant.
  void setClientServerCommunication ( com::PtrCommunication communication );

  com::PtrCommunication getClientServerCommunication() const;

  /// Returns true, if the participant uses a master process.
  bool useMaster();

  void setUseMaster(bool useMaster);

  /**
   * @brief Returns true, if the
   */
  //bool isServer();

  //void setIsServer ( bool value );

private:

  logging::Logger _log{"impl::Participant"};

  static int _participantsSize;

  std::string _name;

  int _id;

  std::vector<PtrWatchPoint> _watchPoints;

  /// Export contexts to export meshes, data, and more.
  std::vector<io::ExportContext> _exportContexts;

  std::vector<action::PtrAction> _actions;

  /// All mesh contexts involved in a simulation, mesh ID == index.
  std::vector<MeshContext*> _meshContexts;

  /// Read mapping contexts used by the participant.
  utils::ptr_vector<MappingContext> _readMappingContexts;

  /// Write mapping contexts used by the participant.
  utils::ptr_vector<MappingContext> _writeMappingContexts;

  /// Mesh contexts used by the participant.
  std::vector<MeshContext*> _usedMeshContexts;

  std::vector<DataContext*> _dataContexts;

  utils::ptr_vector<DataContext> _writeDataContexts;

  utils::ptr_vector<DataContext> _readDataContexts;

  //io::ExportContext _exportContext;

  com::PtrCommunication _clientServerCommunication;

  bool _useMaster = false;

  template<typename ELEMENT_T>
  bool isDataValid (
    const std::vector<ELEMENT_T>& data,
    const ELEMENT_T&              newElement ) const;

  void checkDuplicatedUse ( const mesh::PtrMesh& mesh );

  void checkDuplicatedData ( const mesh::PtrData& data );

  /// To allow white box tests.
  friend struct PreciceTests::Serial::TestConfiguration;
};


// --------------------------------------------------------- HEADER DEFINITIONS


template< typename ELEMENT_T >
bool Participant:: isDataValid
(
  const std::vector<ELEMENT_T>& data,
  const ELEMENT_T&              newElement ) const
{
  for ( size_t i=0; i < data.size(); i++ ) {
    if ( data[i].name == newElement.name ) {
      return false;
    }
  }
  return true;
}

}} // namespace precice, impl
