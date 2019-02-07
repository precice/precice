#include "AxialGeoMultiscaleMapping.hpp"
#include "mesh/Mesh.hpp"

namespace precice {

namespace mapping {

AxialGeoMultiscaleMapping:: AxialGeoMultiscaleMapping
(
  Constraint constraint,
  int        dimensions,
  MultiscaleType type,
  double radius)
:
  Mapping(constraint, dimensions),
  _type(type),
  _radius(radius)
{
  setInputRequirement(Mapping::MeshRequirement::VERTEX);
  setOutputRequirement(Mapping::MeshRequirement::VERTEX);

  assertion(dimensions == 2 || dimensions == 3);
  if(dimensions == 2){
    _scaling = 1.5;
  } else if(dimensions == 3){
    _scaling = 2.0;
  }
}

void AxialGeoMultiscaleMapping:: computeMapping()
{
  TRACE(output()->vertices().size());

  assertion(input().get() != nullptr);
  assertion(output().get() != nullptr);

  if (getConstraint() == CONSISTENT){
    DEBUG("Compute consistent mapping");
    assertion(_type == SPREAD, "Not yet implemented");
    CHECK(input()->vertices().size()==1, "You can only define an axial geometric multiscale mapping of type spread "
        "from a mesh with exactly one vertex. Mesh " << input()->getName() << " has " << input()->vertices().size() << " vertices.");

    // Nothing to do here
  }
  else {
    assertion(getConstraint() == CONSERVATIVE, getConstraint());
    assertion(false, "Not yet implemented");
    DEBUG("Compute conservative mapping");
  }
  _hasComputedMapping = true;
}

bool AxialGeoMultiscaleMapping:: hasComputedMapping() const
{
  TRACE(_hasComputedMapping);
  return _hasComputedMapping;
}

void AxialGeoMultiscaleMapping:: clear()
{
  TRACE();
  _hasComputedMapping = false;
}

void AxialGeoMultiscaleMapping:: map
(
  int inputDataID,
  int outputDataID )
{
  TRACE(inputDataID, outputDataID);

  const Eigen::VectorXd& inputValues = input()->data(inputDataID)->values();
  int valueDimensions = input()->data(inputDataID)->getDimensions();

  assertion(inputValues.size() == valueDimensions);

  Eigen::VectorXd& outputValues = output()->data(outputDataID)->values();

  assertion(input()->vertices().size() == 1);
  mesh::Vertex& v0 = input()->vertices()[0];

  assertion ( valueDimensions == output()->data(outputDataID)->getDimensions(),
              valueDimensions, output()->data(outputDataID)->getDimensions() );
  assertion ( inputValues.size() / valueDimensions == (int)input()->vertices().size(),
               inputValues.size(), valueDimensions, input()->vertices().size() );
  assertion ( outputValues.size() / valueDimensions == (int)output()->vertices().size(),
               outputValues.size(), valueDimensions, output()->vertices().size() );

  if (getConstraint() == CONSISTENT){
    DEBUG("Map consistent");
    assertion(_type == SPREAD, "Not yet implemented");
    size_t const outSize = output()->vertices().size();
    for ( size_t i=0; i < outSize; i++ ){

      Eigen::VectorXd difference(getDimensions());
      difference = v0.getCoords();
      difference -= output()->vertices()[i].getCoords();
      double distance = difference.norm() / _radius;
      CHECK(distance <= 1.05, "Mesh " << output()->getName() << " has vertices that do not coincide with the geometric"
          " multiscale interface define by mesh " << input()->getName() << ". Radius is " << _radius << " and distance between"
          " 2D/3D point and 1D point is " << difference.norm());

      for ( int dim=0; dim < valueDimensions; dim++ ){
        // only parabolic shape supported
        outputValues((i*valueDimensions)+dim) = inputValues(dim) * (1.0 - distance*distance) * _scaling;
      }
    }
  }
  else {
    assertion(getConstraint() == CONSERVATIVE, getConstraint());
    assertion(false, "Not yet implemented");
    DEBUG("Map conservative");
  }
}

void AxialGeoMultiscaleMapping::tagMeshFirstRound()
{
  TRACE();

  computeMapping();

  if (getConstraint() == CONSISTENT){
    assertion(_type == SPREAD, "Not yet implemented");
    assertion( input()->vertices().size() == 1);

    input()->vertices()[0].tag();

  }
  else {
    assertion(getConstraint() == CONSERVATIVE, getConstraint());
    assertion(false, "Not yet implemented");
  }

  clear();
}

void AxialGeoMultiscaleMapping::tagMeshSecondRound()
{
  TRACE();
  // no operation needed here for the moment
}

}} // namespace precice, mapping
