// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IMPL_COUPLINGPARTICIPANT_HPP_
#define PRECICE_IMPL_COUPLINGPARTICIPANT_HPP_

#include "precice/Constants.hpp"
#include "SharedPointer.hpp"
#include "action/SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"
#include "geometry/SharedPointer.hpp"
#include "mapping/SharedPointer.hpp"
#include "spacetree/SharedPointer.hpp"
#include "io/config/ExportConfiguration.hpp"
#include "io/ExportContext.hpp"
#include "com/Communication.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "utils/Helpers.hpp"
#include "tarch/logging/Log.h"
#include "utils/Dimensions.hpp"
#include "utils/PointerVector.hpp"
#include "boost/tuple/tuple.hpp"
#include <map>
#include <string>

namespace precice {
  namespace impl {
    struct DataContext;
    struct MeshContext;
    struct MappingContext;
  }
  namespace com {
    class Communication;
  }
  namespace tests {
    class SolverInterfaceTestGeometry;
    class SolverInterfaceTest;
    class SolverInterfaceTestRemote;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace impl {

/**
 * @brief Holds coupling state of one participating solver in coupled simulation.
 */
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
   * @param name [IN] Name of the participant. Has to be unique.
   */
  Participant (
    const std::string&          name,
    mesh::PtrMeshConfiguration& meshConfig );

  /**
   * @brief Destructor.
   */
  virtual ~Participant();

  /**
   * @brief Returns the name of the participant.
   */
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

  /**
   * @brief Adds a geometry to be used by the participant.
   */
  void useMesh (
    const mesh::PtrMesh&                   mesh,
    const geometry::PtrGeometry&   geometry,
    const spacetree::PtrSpacetree& spacetree,
    const utils::DynVector&                localOffset,
    bool                                   remote,
    const std::string&                     fromParticipant,
    double                                 safetyFactor,
    bool                                   provideMesh );

  void addAction ( const action::PtrAction& action );

  std::vector<action::PtrAction>& actions();

  const std::vector<action::PtrAction>& actions() const;

  /**
   * @brief Adds an export context to export meshes and data.
   */
  void addExportContext ( const io::ExportContext& context );

  /**
   * @brief Returns all export contexts for exporting meshes and data.
   */
  const std::vector<io::ExportContext>& exportContexts() const;

  /**
   * @brief Returns true, if the participant uses a precice in form of a server.
   */
  bool useServer();

  /**
   * @brief Sets the client-server com. for the participant.
   */
  void setClientServerCommunication ( com::Communication::SharedPointer communication );

  com::Communication::SharedPointer getClientServerCommunication() const;

  /**
   * @brief Returns true, if the participant uses a master process.
   */
  bool useMaster();

  /**
   * @brief Sets the masterall com. for the participant.
   */
  void setMasterSlaveCommunication ( com::Communication::SharedPointer communication );

  com::Communication::SharedPointer getMasterSlaveCommunication() const;

  /**
   * @brief Returns true, if the
   */
  //bool isServer();

  //void setIsServer ( bool value );

private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  static int _participantsSize;

  std::string _name;

  int _id;

  std::vector<PtrWatchPoint> _watchPoints;

  // @brief Export contexts to export meshes, data, and more.
  std::vector<io::ExportContext> _exportContexts;

  std::vector<action::PtrAction> _actions;

  // @brief All mesh contexts involved in a simulation, mesh ID == index.
  std::vector<MeshContext*> _meshContexts;

  // @brief Read mapping contexts used by the participant.
  utils::ptr_vector<MappingContext> _readMappingContexts;

  // @brief Write mapping contexts used by the participant.
  utils::ptr_vector<MappingContext> _writeMappingContexts;

  // @brief Mesh contexts used by the participant.
  std::vector<MeshContext*> _usedMeshContexts;

  std::vector<DataContext*> _dataContexts;

  utils::ptr_vector<DataContext> _writeDataContexts;

  utils::ptr_vector<DataContext> _readDataContexts;

  //io::ExportContext _exportContext;


  com::Communication::SharedPointer _clientServerCommunication;

  com::Communication::SharedPointer _masterSlaveCommunication;


  template<typename ELEMENT_T>
  bool isDataValid (
    const std::vector<ELEMENT_T>& data,
    const ELEMENT_T&              newElement ) const;

  void checkDuplicatedUse ( const mesh::PtrMesh& mesh );

  void checkDuplicatedData ( const mesh::PtrData& data );

  // @brief To allow white box tests.
  friend class tests::SolverInterfaceTest;
  friend class tests::SolverInterfaceTestRemote;
  friend class tests::SolverInterfaceTestGeometry;
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

#endif /* PRECICE_IMPL_COUPLINGPARTICIPANT_HPP_ */
