#include "Mapping.hpp"
#include <boost/config.hpp>
#include <ostream>
#include "mesh/Utils.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/assertion.hpp"

namespace precice {
namespace mapping {

Mapping::Mapping(
    Constraint constraint,
    int        dimensions)
    : _constraint(constraint),
      _inputRequirement(MeshRequirement::UNDEFINED),
      _outputRequirement(MeshRequirement::UNDEFINED),
      _input(),
      _output(),
      _dimensions(dimensions)
{
}

void Mapping::setMeshes(
    const mesh::PtrMesh &input,
    const mesh::PtrMesh &output)
{
  _input  = input;
  _output = output;
}

const mesh::PtrMesh &Mapping::getInputMesh() const
{
  return _input;
}

const mesh::PtrMesh &Mapping::getOutputMesh() const
{
  return _output;
}

Mapping::Constraint Mapping::getConstraint() const
{
  return _constraint;
}

Mapping::MeshRequirement Mapping::getInputRequirement() const
{
  return _inputRequirement;
}

Mapping::MeshRequirement Mapping::getOutputRequirement() const
{
  return _outputRequirement;
}

mesh::PtrMesh Mapping::input() const
{
  return _input;
}

mesh::PtrMesh Mapping::output() const
{
  return _output;
}

void Mapping::setInputRequirement(
    MeshRequirement requirement)
{
  _inputRequirement = requirement;
}

void Mapping::setOutputRequirement(
    MeshRequirement requirement)
{
  _outputRequirement = requirement;
}

int Mapping::getDimensions() const
{
  return _dimensions;
}

void Mapping::scaleConsistentMapping(int inputDataID, int outputDataID) const
{
  // Only serial participant is supported for scale-consistent mapping
  PRECICE_ASSERT((not utils::MasterSlave::isMaster()) and (not utils::MasterSlave::isSlave()));

  // If rank is not empty and do not contain connectivity information, raise error
  if ((input()->edges().empty() and (not input()->vertices().empty())) or
      (((input()->getDimensions() == 3) and input()->triangles().empty()) and (not input()->vertices().empty()))) {
    logging::Logger _log{"mapping::Mapping"};
    PRECICE_ERROR("Connectivity information is missing for the mesh {}. "
                  "Scaled consistent mapping requires connectivity information.",
                  input()->getName());
  }
  if ((output()->edges().empty() and (not output()->vertices().empty())) or
      (((output()->getDimensions() == 3) and output()->triangles().empty()) and (not output()->vertices().empty()))) {
    logging::Logger _log{"mapping::Mapping"};
    PRECICE_ERROR("Connectivity information is missing for the mesh {}. "
                  "Scaled consistent mapping requires connectivity information.",
                  output()->getName());
  }

  auto &outputValues    = output()->data(outputDataID)->values();
  int   valueDimensions = input()->data(inputDataID)->getDimensions();

  // Integral is calculated on each direction separately
  auto integralInput  = mesh::integrate(input(), input()->data(inputDataID));
  auto integralOutput = mesh::integrate(output(), output()->data(outputDataID));

  // Create reshape the output values vector to matrix
  Eigen::Map<Eigen::MatrixXd> outputValuesMatrix(outputValues.data(), valueDimensions, outputValues.size() / valueDimensions);

  // Scale in each direction
  Eigen::VectorXd scalingFactor = integralInput.array() / integralOutput.array();
  outputValuesMatrix.array().colwise() *= scalingFactor.array();
}

bool Mapping::hasConstraint(const Constraint &constraint) const
{
  return (getConstraint() == constraint);
}

bool operator<(Mapping::MeshRequirement lhs, Mapping::MeshRequirement rhs)
{
  switch (lhs) {
  case (Mapping::MeshRequirement::UNDEFINED):
    return rhs != Mapping::MeshRequirement::UNDEFINED;
  case (Mapping::MeshRequirement::VERTEX):
    return rhs == Mapping::MeshRequirement::FULL;
  case (Mapping::MeshRequirement::FULL):
    return false;
  };
  BOOST_UNREACHABLE_RETURN(false);
}

std::ostream &operator<<(std::ostream &out, Mapping::MeshRequirement val)
{
  switch (val) {
  case (Mapping::MeshRequirement::UNDEFINED):
    out << "UNDEFINED";
    break;
  case (Mapping::MeshRequirement::VERTEX):
    out << "VERTEX";
    break;
  case (Mapping::MeshRequirement::FULL):
    out << "FULL";
    break;
  default:
    PRECICE_ASSERT(false, "Implementation does not cover all cases");
  };
  return out;
}

} // namespace mapping
} // namespace precice
