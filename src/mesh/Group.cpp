#include "Group.hpp"

namespace precice {
namespace mesh {

void Group:: add
(
  Vertex& vertex )
{
  _vertices.push_back ( &vertex );
}

void Group:: add
(
  Vertex* vertex )
{
  PRECICE_ASSERT( vertex != nullptr );
  _vertices.push_back ( vertex );
}

void Group:: add
(
  Edge& edge )
{
  _edges.push_back ( &edge );
}

void Group:: add
(
  Edge* edge )
{
  PRECICE_ASSERT(edge != nullptr);
  _edges.push_back(edge);
}

void Group:: add
(
  Triangle& triangle )
{
  _triangles.push_back(&triangle);
}
void Group:: add
(
  Triangle* triangle )
{
  PRECICE_ASSERT( triangle != nullptr );
  _triangles.push_back ( triangle );
}

void Group:: add
(
  Quad& quad )
{
  _quads.push_back(&quad);
}
void Group:: add
(
  Quad* quad )
{
  PRECICE_ASSERT( quad != nullptr );
  _quads.push_back ( quad );
}

void Group:: add
(
  Group& group )
{
  _vertices.insert ( _vertices.end(), group.vertices().begin(),
                     group.vertices().end() );
  _edges.insert ( _edges.end(), group.edges().begin(), group.edges().end() );
  _triangles.insert ( _triangles.end(), group.triangles().begin(),
                      group.triangles().end() );
}

Group::VertexContainer& Group:: vertices()
{
  return _vertices;
}

const Group::VertexContainer& Group:: vertices() const
{
  return _vertices;
}

Group::EdgeContainer& Group:: edges()
{
  return _edges;
}

const Group::EdgeContainer& Group:: edges() const
{
  return _edges;
}

Group::TriangleContainer& Group:: triangles()
{
  return _triangles;
}

const Group::TriangleContainer& Group:: triangles() const
{
  return _triangles;
}

Group::QuadContainer& Group:: quads()
{
  return _quads;
}

const Group::QuadContainer& Group:: quads() const
{
  return _quads;
}

bool Group:: empty() const
{
  return _vertices.empty() && _edges.empty() && _triangles.empty() && _quads.empty();
}

size_t Group:: size() const
{
  return _vertices.size() + _edges.size() + _triangles.size() + _quads.size();
}

void Group:: clear()
{
  _vertices.clear();
  _edges.clear();
  _triangles.clear();
  _quads.clear();
}

}} // namespace precice, mesh
