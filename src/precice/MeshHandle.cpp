// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "precice/MeshHandle.hpp"
#include "utils/Dimensions.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Group.hpp"

namespace precice {
  namespace impl {

    struct VertexIteratorImplementation {
      mesh::Group::VertexContainer::const_iterator iterator;
      //double coords[];
      //int dimensions;
    };

    struct EdgeIteratorImplementation {
      mesh::Group::EdgeContainer::const_iterator iterator;
      //double coords[];
      //int dimensions;
    };

    struct TriangleIteratorImplementation {
      mesh::Group::TriangleContainer::const_iterator iterator;
      //double coords[];
      //int dimensions;
    };
  }
}

namespace precice {

VertexIterator:: VertexIterator
(
  const mesh::Group& content,
  bool               begin )
:
  _impl (new impl::VertexIteratorImplementation())
{
  if ( begin ) {
    _impl->iterator = content.vertices().begin();
  }
  else  {
    _impl->iterator = content.vertices().end();
  }
  //_impl->dimensions = content.getDimensions();
}

VertexIterator:: VertexIterator
(
  const VertexIterator& toCopy )
:
  _impl ( new impl::VertexIteratorImplementation(*toCopy._impl) )
{}

VertexIterator & VertexIterator:: operator=
(
  const VertexIterator& toAssign )
{
  _impl->iterator = toAssign._impl->iterator;
//  assertion2 ( _impl->dimensions == toAssign._impl->dimensions,
//               _impl->dimensions, toAssign._impl->dimensions );
//  for ( int dim=0; dim < _impl->dimensions; dim++ ){
//    _impl->coords[dim] = toAssign._impl->coords[dim];
//  }
  return *this;
}

VertexIterator:: ~VertexIterator()
{
  assertion ( _impl != NULL );
  delete _impl;
}

VertexIterator& VertexIterator:: operator++(int)
{
  _impl->iterator++;
  return *this;
}


const double* VertexIterator:: vertexCoords()
{
//  for ( int dim=0; dim < _impl->dimensions; dim++ ){
//    _impl->coords[dim] = (*_impl->iterator).getCoords()[dim];
//  }
//  return _impl->coords;
  return tarch::la::raw((*_impl->iterator).getCoords());
}

int VertexIterator:: vertexID()
{
  return (*_impl->iterator).getID();
}

bool VertexIterator:: operator!=
(
  const VertexIterator& vertexIterator )
{
  return _impl->iterator != vertexIterator._impl->iterator;
}


VertexHandle:: VertexHandle
(
  const mesh::Group& content )
:
  _content ( content )
{}

VertexIterator VertexHandle:: begin() const
{
  return VertexIterator(_content, true);
}

VertexIterator VertexHandle:: end() const
{
  return VertexIterator(_content, false);
}

size_t VertexHandle:: size() const
{
  return _content.vertices().size();
}


EdgeIterator:: EdgeIterator
(
  const mesh::Group& content,
  bool               begin )
:
  _impl (new impl::EdgeIteratorImplementation())
{
  if ( begin ) {
    _impl->iterator = content.edges().begin();
  }
  else  {
    _impl->iterator = content.edges().end();
  }
}

EdgeIterator:: ~EdgeIterator()
{
  assertion ( _impl != NULL );
  delete _impl;
}

EdgeIterator& EdgeIterator:: operator++(int)
{
  _impl->iterator++;
  return *this;
}

const double* EdgeIterator:: vertexCoords
(
  int vertexIndex )
{
  //using tarch::la::dwrap;
  //dwrap(_impl->coords) = (*_impl->iterator).vertex(vertexIndex).getCoords();
  const utils::DynVector& coords = (*_impl->iterator).vertex(vertexIndex).getCoords();
  return tarch::la::raw(coords);
}

int EdgeIterator:: vertexID
(
  int vertexIndex )
{
  return (*_impl->iterator).vertex(vertexIndex).getID();
}

bool EdgeIterator:: operator!=
(
  const EdgeIterator & edgeIterator )
{
  return _impl->iterator != edgeIterator._impl->iterator;
}


EdgeHandle:: EdgeHandle
(
  const mesh::Group & content )
:
  _content ( content )
{}

EdgeIterator EdgeHandle:: begin () const
{
  return EdgeIterator ( _content, true );
}

EdgeIterator EdgeHandle:: end () const
{
  return EdgeIterator ( _content, false );
}

size_t EdgeHandle:: size () const
{
  return _content.edges().size();
}


TriangleIterator:: TriangleIterator
(
  const mesh::Group & content,
  bool                begin )
:
  _impl (new impl::TriangleIteratorImplementation())
{
  if ( begin ) {
    _impl->iterator = content.triangles().begin();
  }
  else  {
    _impl->iterator = content.triangles().end();
  }
}

TriangleIterator:: ~TriangleIterator ()
{
  assertion ( _impl != NULL );
  delete _impl;
}

TriangleIterator & TriangleIterator:: operator++ (int)
{
  _impl->iterator++;
  return *this;
}

const double* TriangleIterator:: vertexCoords
(
  int vertexIndex )
{
//  using tarch::la::dwrap;
//  dwrap(_impl->coords) = (*_impl->iterator).vertex(vertexIndex).getCoords();
//  return _impl->coords;
  const utils::DynVector& coords = (*_impl->iterator).vertex(vertexIndex).getCoords();
  return tarch::la::raw(coords);
}

int TriangleIterator:: vertexID
(
  int vertexIndex )
{
  return (*_impl->iterator).vertex(vertexIndex).getID();
}

bool TriangleIterator:: operator!=
(
  const TriangleIterator & triangleIterator )
{
  return _impl->iterator != triangleIterator._impl->iterator;
}

TriangleHandle:: TriangleHandle
(
  const mesh::Group & content )
:
  _content ( content )
{}


TriangleIterator TriangleHandle:: begin () const
{
  return TriangleIterator ( _content, true );
}


TriangleIterator TriangleHandle:: end () const
{
  return TriangleIterator ( _content, false );
}

size_t TriangleHandle:: size () const
{
  return _content.triangles().size();
}

MeshHandle:: MeshHandle
(
  const mesh::Group & content )
:
  _vertexHandle ( content ),
  _edgeHandle ( content ),
  _triangleHandle ( content )
{}

const VertexHandle & MeshHandle:: vertices () const
{
  return _vertexHandle;
}

const EdgeHandle & MeshHandle:: edges () const
{
  return _edgeHandle;
}

const TriangleHandle & MeshHandle:: triangles () const
{
  return _triangleHandle;
}

} // nammespace precice
