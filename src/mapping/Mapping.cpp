#include "Mapping.hpp"
#include <boost/config.hpp>
#include <ostream>
#include "utils/MasterSlave.hpp"
#include "utils/assertion.hpp"
#include "mesh/Utils.hpp"

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
  // If rank is not empty and do not contain connectivity information, raise error
  if (input()->edges().empty() and (not input()->vertices().empty())) {
    logging::Logger _log{"mapping::Mapping"};
    PRECICE_ERROR("There is no connectivity information defined for mesh " << input()->getName() << ". Scaled consistent mapping requires connectivity information.");
  }
  if (output()->edges().empty() and (not output()->vertices().empty())) {
    logging::Logger _log{"mapping::Mapping"};
    PRECICE_ERROR("There is no connectivity information defined for mesh " << output()->getName() << ". Scaled consistent mapping requires connectivity information.");
  }

  const auto &inputValues  = input()->data(inputDataID)->values();
  auto &      outputValues = output()->data(outputDataID)->values();

  int valueDimensions = input()->data(inputDataID)->getDimensions();
  int meshDimensions  = input()->getDimensions();

  // Integral is calculated on each direction separately
  std::vector<double> integralInput = mesh::integrateOverlap(input(), input()->data(inputDataID));
  std::vector<double> integralOutput = mesh::integrate(output(), output()->data(outputDataID));
  
  // If the mesh is distributed, we need to calculate the global integral
  if (utils::MasterSlave::isMaster() or utils::MasterSlave::isSlave()) {
    std::vector<double> globalInputIntegral(valueDimensions);
    std::vector<double> globalOutputIntegral(valueDimensions);
    utils::MasterSlave::allreduceSum(integralInput.data(), globalInputIntegral.data(), integralInput.size());
    utils::MasterSlave::allreduceSum(integralOutput.data(), globalOutputIntegral.data(), integralOutput.size());
    integralInput  = globalInputIntegral;
    integralOutput = globalOutputIntegral;
  }

  // Scale in each dimension
  size_t const outSize = output()->vertices().size();
  for (int dim = 0; dim < valueDimensions; ++dim) {
    PRECICE_ASSERT(std::fabs(integralOutput.at(dim)) > 1e-9, "Surface integral is calculated as zero in consistent mapping scaling.");
    PRECICE_ASSERT(std::fabs(integralInput.at(dim)) > 1e-9, "Surface integral is calculated as zero in consistent mapping scaling.");
    double scalingFactor = integralInput.at(dim) / integralOutput.at(dim);
    // Scaling factor is different in each direction, cannot directly multiply
    for (size_t i = 0; i < outSize; i++) {
      outputValues((i * valueDimensions) + dim) *= scalingFactor;
    }
  }
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
