#include "Vertex.hpp"
#include "utils/ManageUniqueIDs.hpp"

namespace precice {
namespace mesh {

//Vertex:: Vertex
//(
//  const utils::Vector3D& coordinates,
//  int                    id )
//:
//  PropertyContainer (),
//  _id ( id ),
//  _coords ( coordinates ),
//  _normal ( 3, 0.0 ),
//  _mesh ( NULL )
//{}
//
//Vertex:: Vertex
//(
//  const utils::Vector2D& coordinates,
//  int                    id )
//:
//  PropertyContainer (),
//  _id ( id ),
//  _coords ( coordinates ),
//  _normal ( 2, 0.0 ),
//  _mesh ( NULL )
//{}
//
//Vertex:: Vertex
//(
//  const utils::DynVector& coordinates,
//  int                     id )
//:
//  PropertyContainer (),
//  _id ( id ),
//  _coords ( coordinates ),
//  _normal ( _coords.size(), 0.0 ),
//  _mesh ( NULL )
//{}

//Vertex:: Vertex
//(
//  const utils::Vector3D& coordinates,
//  int                    id,
//  Mesh&                  mesh )
//:
//  PropertyContainer (),
//  _id ( id ),
//  _coords ( coordinates ),
//  _normal ( 3, 0.0 ),
//  _mesh ( & mesh )
//{}
//
//Vertex:: Vertex
//(
//  const utils::Vector2D& coordinates,
//  int                    id,
//  Mesh&                  mesh )
//:
//  PropertyContainer (),
//  _id ( id ),
//  _coords ( coordinates ),
//  _normal ( 2, 0.0 ),
//  _mesh ( & mesh )
//{}
//
//Vertex:: Vertex
//(
//  const utils::DynVector& coordinates,
//  int                     id,
//  Mesh&                   mesh )
//:
//  PropertyContainer (),
//  _id ( id ),
//  _coords ( coordinates ),
//  _normal ( _coords.size(), 0.0 ),
//  _mesh ( & mesh )
//{}

int Vertex:: getDimensions() const
{
  return _coords.size();
}

//void Vertex:: setCoords
//(
//  const utils::DynVector& coordinates )
//{
//  assertion ( coordinates.size() == _coords.size(), coordinates.size(), _coords.size() );
//  _coords = coordinates;
//}
//
//void Vertex:: setNormal
//(
//  const utils::DynVector& normal )
//{
//  assertion ( normal.size() == _normal.size(), normal.size(), _normal.size() );
//  _normal = normal;
//}

const Eigen::VectorXd& Vertex::getNormal () const
{
  return _normal;
}

const Mesh* Vertex:: mesh () const
{
  return _mesh;
}

Mesh* Vertex:: mesh ()
{
  return _mesh;
}

int Vertex:: getGlobalIndex() const {
  return _globalIndex;
}

void Vertex:: setGlobalIndex(int globalIndex){
  _globalIndex = globalIndex;
}

bool Vertex:: isOwner() const {
  return _owner;
}

void Vertex:: setOwner(bool owner){
  _owner = owner;
}

}} // namespace precice, mesh
