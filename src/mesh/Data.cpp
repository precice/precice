// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "Data.hpp"
#include "Mesh.hpp"
#include "utils/Globals.hpp"

namespace precice {
namespace mesh {

tarch::logging::Log Data:: _log ( "precice::mesh::Data" );

size_t Data:: _dataCount = 0;

Data:: Data()
:
  _name ( "" ),
  _id ( -1 ),
  _dimensions ( 0 )
{
  assertion ( false );
}

Data:: Data
(
  const std::string& name,
  int                id,
  int                dimensions )
:
  _values(),
  _name ( name ),
  _id ( id ),
  _dimensions ( dimensions )
{
  assertion1 ( dimensions > 0, dimensions );
  _dataCount ++;
}

Data:: ~Data()
{
  _dataCount --;
}

utils::DynVector& Data:: values()
{
  return _values;
}

const utils::DynVector& Data:: values() const
{
  return _values;
}

const std::string& Data:: getName() const
{
  return _name;
}

int Data:: getID() const
{
  return _id;
}

int Data:: getDimensions() const
{
   return _dimensions;
}

//void Data:: setMesh
//(
//  Mesh* mesh )
//{
//  assertion ( mesh != NULL );
//  _mesh = mesh;
//}

Mesh* Data:: mesh()
{
  return _mesh;
}

const Mesh * Data:: mesh() const
{
  return _mesh;
}

size_t Data:: getDataCount()
{
  return _dataCount;
}

void Data:: resetDataCount()
{
  _dataCount = 0;
}

}} // namespace precice, mesh
