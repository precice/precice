#include "Data.hpp"
#include "Mesh.hpp"

namespace precice {
namespace mesh {

size_t Data:: _dataCount = 0;

Data:: Data()
:
  _name ( "" ),
  _id ( -1 ),
  _dimensions ( 0 )
{
  P_assertion( false );
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
  P_assertion( dimensions > 0, dimensions );
  _dataCount ++;
}

Data:: ~Data()
{
  _dataCount --;
}

Eigen::VectorXd& Data:: values()
{
  return _values;
}

const Eigen::VectorXd& Data:: values() const
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

size_t Data:: getDataCount()
{
  return _dataCount;
}

void Data:: resetDataCount()
{
  _dataCount = 0;
}

}} // namespace precice, mesh
