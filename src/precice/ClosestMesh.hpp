// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CLOSESTMESH_HPP_
#define PRECICE_CLOSESTMESH_HPP_

#include <vector>
#include "Constants.hpp"

namespace precice {
  namespace impl {
    class ClosestMeshImplementation;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {

/**
 * @brief Holds information about a query point in relation to geometry.
 *
 * Is returned from the function inquireClosestMesh().
 */
struct ClosestMesh
{
public:

  /**
   * @brief Constructor.
   */
  ClosestMesh ( int dimensions );

  /**
   * @brief Copy constructor.
   */
  ClosestMesh ( const ClosestMesh& toCopy );

  /**
   * @brief Destructor.
   */
  ~ClosestMesh();

  /**
   * @brief Assignment operator.
   */
  ClosestMesh& operator= ( const ClosestMesh& toAssign );

  /**
   * @brief Returns the mesh IDs assigned to the closest mesh.
   */
  std::vector<int> & meshIDs();

  /**
   * @brief Returns the relative position of the closest mesh.
   *
   * The meaning of the position can be checked by using the
   * precice::constants::position...() methods.
   */
  int position();

  /**
   * @brief Sets the position.
   */
  void setPosition ( int position );

  /**
   * @brief Returns the distance vector to the closest mesh.
   *
   * This vector will (in suitable cases) be orthogonal to an edge or triangle
   * of the closest mesh. However, in cases where it is not possible to find an
   * orthogonal projection onto an edge/triangle, the vector will point to the
   * closest vertex (which is an approximation to the normal in some sense).
   */
  const double* distanceVector ();

  /**
   * @brief Sets the distance vector (and distance).
   */
  void setDistanceVector ( const double* distanceVector );

  /**
   * @brief Returns the distance (unsigned) to the closest mesh.
   */
  double distance();

private:

  impl::ClosestMeshImplementation* _impl;
};

} // namespace precice

#endif /* PRECICE_CLOSESTMESH_HPP_ */
