// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "FindClosestVertex.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Mesh.hpp"

namespace precice {
namespace query {

const utils::DynVector& FindClosestVertex:: getSearchPoint() const
{
  return _searchPoint;
}

bool FindClosestVertex:: hasFound() const
{
  return _closestVertex != NULL;
}

double FindClosestVertex:: getEuclidianDistance()
{
   return _shortestDistance;
}

mesh::Vertex& FindClosestVertex:: getClosestVertex()
{
  assertion ( _closestVertex != NULL );
  return *_closestVertex;
}

}} // namespace precice, query
