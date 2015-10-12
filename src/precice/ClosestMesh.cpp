// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ClosestMesh.hpp"
#include "utils/Dimensions.hpp"
#include "spacetree/Spacetree.hpp"
#include <limits>

namespace precice {
  namespace impl {

    /**
     * @brief Holds data of ClosestMesh to hide it from interface of ClosestMesh.
     */
    struct ClosestMeshImplementation
    {
      // @brief Geometry IDs of the closest geometry.
      std::vector<int> meshIDs;

      // @brief Position of the query point relative to the closest geometry.
      int position;

      // @brief Shortest distance vector to closest mesh.
      utils::DynVector distanceVector;

      /**
       * @brief Constructor.
       */
      ClosestMeshImplementation ( int dimensions )
      :
        meshIDs(),
        position ( spacetree::Spacetree::positionUndefined() ),
        distanceVector ( dimensions )
      {}
    };

  }
}

namespace precice {

ClosestMesh:: ClosestMesh
(
  int dimensions )
:
  _impl ( new impl::ClosestMeshImplementation(dimensions) )
{

  _impl->position = spacetree::Spacetree::positionUndefined();
  assign(_impl->distanceVector) = std::numeric_limits<double>::max();
}

ClosestMesh:: ClosestMesh ( const ClosestMesh& toCopy )
:
  _impl ( new impl::ClosestMeshImplementation(toCopy._impl->distanceVector.size()) )
{
  _impl->meshIDs = toCopy._impl->meshIDs;
  _impl->position = toCopy._impl->position;
  _impl->distanceVector = toCopy._impl->distanceVector;
}

ClosestMesh:: ~ClosestMesh()
{
  assertion ( _impl != nullptr );
  delete _impl;
}

ClosestMesh& ClosestMesh:: operator=
(
  const ClosestMesh& toAssign )
{
  _impl->meshIDs = toAssign._impl->meshIDs;
  _impl->position = toAssign._impl->position;
  _impl->distanceVector = toAssign._impl->distanceVector;
  return *this;
}

std::vector<int>& ClosestMesh:: meshIDs()
{
  assertion ( _impl != nullptr );
  return _impl->meshIDs;
}

int ClosestMesh:: position()
{
  assertion ( _impl != nullptr );
  return _impl->position;
}

void ClosestMesh:: setPosition
(
  int position )
{
  assertion ( _impl != nullptr );
  _impl->position = position;
}

const double* ClosestMesh:: distanceVector()
{
  assertion ( _impl != nullptr );
  return tarch::la::raw(_impl->distanceVector);
}

void ClosestMesh:: setDistanceVector
(
  const double* distanceVector )
{
  assertion ( _impl != nullptr );
  assertion ( distanceVector != nullptr );
  const int dim = _impl->distanceVector.size();
  for ( int i=0; i < dim; i++ ){
    _impl->distanceVector[i] = distanceVector[i];
  }
}

double ClosestMesh:: distance ()
{
  assertion ( _impl != nullptr );
  return tarch::la::norm2(_impl->distanceVector);
}

} // namespace precice
