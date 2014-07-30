// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_IMPL_SHAREDPOINTER_HPP_
#define PRECICE_IMPL_SHAREDPOINTER_HPP_

#include "boost/smart_ptr.hpp"

namespace precice {
namespace impl {

class Participant;
typedef boost::shared_ptr<Participant> PtrParticipant;

class Coupling;
typedef boost::shared_ptr<Coupling> PtrCoupling;

class WatchPoint;
typedef boost::shared_ptr<WatchPoint> PtrWatchPoint;

class AbstractDataAction;
typedef boost::shared_ptr<AbstractDataAction> PtrAbstractDataAction;

struct MeshContext;
typedef boost::shared_ptr<MeshContext> PtrMeshContext;

// TODO delete if really not requierd
//class AbstractRequest;
//typedef boost::shared_ptr<AbstractRequest> PtrAbstractRequest;

//class MappingRequest;
//typedef boost::shared_ptr<MappingRequest> PtrMappingRequest;

//class MapDoubleToCouplingMeshRequest;
//typedef boost::shared_ptr<MapDoubleToCouplingMeshRequest>
//                                                   PtrMapToCouplingMeshRequest;

//class Client;
//typedef boost::shared_ptr<Client> PtrClient;
//
//class Server;
//typedef boost::shared_ptr<Server> PtrServer;



}} // namespace precice, impl


#endif /* PRECICE_IMPL_SHAREDPOINTER_HPP_ */
