// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_COM_COMMUNICATEMESH_HPP_
#define PRECICE_COM_COMMUNICATEMESH_HPP_

#include "tarch/logging/Log.h"
#include "mesh/Mesh.hpp"

namespace precice {
   namespace mesh {
      class Mesh;
   }
}

// ------------------------------------------------------------ CLASS DEFINTION

namespace precice {
namespace com {

/**
 * @brief Copies a Mesh object from a sender to a receiver.
 */
class CommunicateMesh
{
public:

  /**
   * @brief Constructor, takes communication to be used in transfer.
   */
  CommunicateMesh (
    com::Communication::SharedPointer communication );

  /**
   * @brief Sends a constructed CustomGeometry to the receiver with given rank.
   */
  void sendMesh (
    const mesh::Mesh & mesh,
    int                rankReceiver );

  /**
   * @brief Copies a CustomGeometry from the sender with given rank.
   */
  void receiveMesh (
    mesh::Mesh & mesh,
    int          rankSender );

  void sendBoundingBox (
    const mesh::Mesh::BoundingBox & bb,
    int                rankReceiver );

  void receiveBoundingBox (
    mesh::Mesh::BoundingBox & bb,
    int          rankSender );


private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  // @brief Communication means used for the transfer of the geometry.
  com::Communication::SharedPointer _communication;
};

}} // namespace precice, com

#endif /* PRECICE_COM_COMMUNICATEMESH_HPP_ */
