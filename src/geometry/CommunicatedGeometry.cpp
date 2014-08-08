// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "CommunicatedGeometry.hpp"
#include "com/CommunicateMesh.hpp"
#include "com/Communication.hpp"
#include "mesh/Mesh.hpp"
#include "utils/Globals.hpp"
#include "utils/Helpers.hpp"

namespace precice {
namespace geometry {

tarch::logging::Log CommunicatedGeometry:: _log ( "precice::geometry::CommunicatedGeometry" );

CommunicatedGeometry:: CommunicatedGeometry
(
  const utils::DynVector&  offset,
  const std::string&       accessor,
  const std::string&       provider )
:
  Geometry ( offset ),
  _accessorName ( accessor ),
  _providerName ( provider ),
  _receivers ()
{
  preciceTrace2 ( "CommunicatedGeometry()", accessor, provider );
}

void CommunicatedGeometry:: addReceiver
(
  const std::string&     receiver,
  com::PtrCommunication com )
{
  preciceTrace1 ( "addReceiver()", receiver );
  assertion ( com.get() != NULL );
  preciceCheck ( ! utils::contained(receiver, _receivers),
                 "addReceiver()", "Receiver \"" << receiver
                 << "\" has been added already to communicated geometry!" );
  preciceCheck ( receiver != _providerName, "addReceiver()",
                 "Receiver \"" << receiver << "\" cannot be the same as "
                 << "provider in communicated geometry!" );
  _receivers[receiver] = com;
}

void CommunicatedGeometry:: specializedCreate
(
  mesh::Mesh& seed )
{
  preciceTrace1 ( "specializedCreate()", seed.getName() );
  preciceCheck ( _receivers.size() > 0, "specializedCreate()",
                 "No receivers specified for communicated geometry to create "
                 << "mesh \"" << seed.getName() << "\"!" );
  if ( _accessorName == _providerName ) {
    preciceCheck ( seed.vertices().size() > 0,
                   "specializedCreate()", "Participant \"" << _accessorName
                   << "\" provides an invalid (possibly empty) mesh \""
                   << seed.getName() << "\"!" );
    typedef std::map<std::string,com::PtrCommunication>::value_type Pair;
    foreach ( Pair & pair, _receivers ) {
      if ( ! pair.second->isConnected() ) {
        pair.second->acceptConnection ( _providerName, pair.first, 0, 1 );
      }
      com::CommunicateMesh(pair.second).sendMesh ( seed, 0 );
    }
  }
  else if ( utils::contained(_accessorName, _receivers) ) {
    assertion ( seed.vertices().size() == 0 );
    assertion ( utils::contained(_accessorName, _receivers) );
    com::PtrCommunication com ( _receivers[_accessorName] );
    if ( ! com->isConnected() ) {
      com->requestConnection ( _providerName, _accessorName, 0, 1 );
    }
    com::CommunicateMesh(com).receiveMesh ( seed, 0 );
  }
  else {
    preciceError( "specializedCreate()", "Participant \"" << _accessorName
                  << "\" uses a communicated geometry to create mesh \""
                  << seed.getName()
                  << "\" but is neither provider nor receiver!" );
  }
}

}} // namespace precice, geometry
