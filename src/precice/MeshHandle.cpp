#include "precice/MeshHandle.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Group.hpp"

namespace precice {
  namespace impl {

    struct VertexIteratorImplementation {
      mesh::Group::VertexContainer::const_iterator iterator;
    };

    struct EdgeIteratorImplementation {
      mesh::Group::EdgeContainer::const_iterator iterator;
    };

    struct TriangleIteratorImplementation {
      mesh::Group::TriangleContainer::const_iterator iterator;
    };
  }
}

namespace precice {

// VertexIterator

VertexIterator::VertexIterator()
    : _impl(NULL) {}

VertexIterator::~VertexIterator()
{
  delete _impl;
}

VertexIterator:: VertexIterator
(
  const mesh::Group& content,
  bool               begin )
{
  if ( begin ) {
    _impl = new Impl{content.vertices().begin()};
  }
  else  {
    _impl = new Impl{content.vertices().end()};
  }
}

VertexIterator:: VertexIterator
(
  const VertexIterator& toCopy ) : _impl(NULL)
{
    if (toCopy._impl) {
        this->_impl = new Impl{toCopy._impl->iterator};
    }
}

VertexIterator & VertexIterator:: operator=
(
  const VertexIterator& toAssign )
{
    VertexIterator cpy(toAssign);
    using std::swap;
    swap(*this, cpy);
    return *this;
}

VertexIterator VertexIterator:: operator++(int unused)
{
  VertexIterator cpy(*this);
  this->operator++();
  return cpy;
}

VertexIterator& VertexIterator::operator++()
{
  _impl->iterator.operator++();
  return *this;
}

const VertexIterator VertexIterator::operator*() const
{
  return *this;
}

bool VertexIterator::operator==(const VertexIterator& other) const
{
  return _impl->iterator == other._impl->iterator;
}

bool VertexIterator::operator!=(const VertexIterator& other) const
{
  return !((*this) == other);
}

const double* VertexIterator:: vertexCoords() const
{
  return (*_impl->iterator).getCoords().data();
}

int VertexIterator:: vertexID() const
{
  return (*_impl->iterator).getID();
}

void VertexIterator::swap(VertexIterator& other)
{
    std::swap(_impl, other._impl);
}

void swap(VertexIterator& lhs, VertexIterator& rhs)
{
    lhs.swap(rhs);
}

// VertexHandle

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

std::size_t VertexHandle:: size() const
{
  return _content.vertices().size();
}

// EdgeIterator

EdgeIterator::EdgeIterator()
    : _impl(NULL) {}

EdgeIterator::~EdgeIterator()
{
  delete _impl;
}

EdgeIterator:: EdgeIterator
(
  const mesh::Group& content,
  bool               begin )
{
  if ( begin ) {
    _impl=new Impl{content.edges().begin()};
  }
  else  {
    _impl=new Impl{content.edges().end()};
  }
}

EdgeIterator:: EdgeIterator(const EdgeIterator& toCopy) : _impl(NULL)
{
    if (toCopy._impl) {
        this->_impl = new Impl{toCopy._impl->iterator};
    }
}

EdgeIterator & EdgeIterator:: operator=(const EdgeIterator& other)
{
    EdgeIterator cpy(other);
    using std::swap;
    swap(*this, cpy);
    return *this;
}

EdgeIterator EdgeIterator:: operator++(int unused)
{
  EdgeIterator cpy(*this);
  this->operator++();
  return cpy;
}

EdgeIterator& EdgeIterator::operator++()
{
  _impl->iterator.operator++();
  return *this;
}

const EdgeIterator EdgeIterator::operator*() const
{
  return *this;
}

bool EdgeIterator::operator==(const EdgeIterator& other) const
{
  return _impl->iterator == other._impl->iterator;
}

bool EdgeIterator::operator!=(const EdgeIterator& other) const
{
  return !((*this) == other);
}


const double* EdgeIterator:: vertexCoords(int vertexIndex) const
{
  return (*_impl->iterator).vertex(vertexIndex).getCoords().data();
}

int EdgeIterator:: vertexID(int vertexIndex) const
{
  return (*_impl->iterator).vertex(vertexIndex).getID();
}

void EdgeIterator::swap(EdgeIterator& other)
{
    std::swap(_impl, other._impl);
}

void swap(EdgeIterator& lhs, EdgeIterator& rhs)
{
    lhs.swap(rhs);
}

// EdgeHandle

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

std::size_t EdgeHandle:: size () const
{
  return _content.edges().size();
}

// TriangleIterator

TriangleIterator::TriangleIterator()
    : _impl(NULL) {}

TriangleIterator::~TriangleIterator()
{
  delete _impl;
}

TriangleIterator:: TriangleIterator
(
  const mesh::Group& content,
  bool               begin )
{
  if ( begin ) {
    _impl = new Impl{content.triangles().begin()};
  }
  else  {
    _impl = new Impl{content.triangles().end()};
  }
}

TriangleIterator:: TriangleIterator(const TriangleIterator& toCopy) : _impl(NULL)
{
    if (toCopy._impl) {
        this->_impl = new Impl{toCopy._impl->iterator};
    }
}

TriangleIterator & TriangleIterator:: operator=(const TriangleIterator& other)
{
    TriangleIterator cpy(other);
    using std::swap;
    swap(*this, cpy);
    return *this;
}

TriangleIterator TriangleIterator:: operator++(int unused)
{
  TriangleIterator cpy(*this);
  this->operator++();
  return cpy;
}

TriangleIterator& TriangleIterator::operator++()
{
  _impl->iterator.operator++();
  return *this;
}

const TriangleIterator TriangleIterator::operator*() const
{
  return *this;
}

bool TriangleIterator::operator==(const TriangleIterator& other) const
{
  return _impl->iterator == other._impl->iterator;
}

bool TriangleIterator::operator!= ( const TriangleIterator& other ) const
{
  return !((*this) == other);
}


const double* TriangleIterator:: vertexCoords ( int vertexIndex ) const
{
  return (*_impl->iterator).vertex(vertexIndex).getCoords().data();
}

int TriangleIterator:: vertexID ( int vertexIndex ) const
{
  return (*_impl->iterator).vertex(vertexIndex).getID();
}

void TriangleIterator::swap(TriangleIterator& other)
{
    std::swap(_impl, other._impl);
}

void swap(TriangleIterator& lhs, TriangleIterator& rhs)
{
    lhs.swap(rhs);
}

// TriangleHandle

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

std::size_t TriangleHandle:: size () const
{
  return _content.triangles().size();
}

// MeshHandle

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
