#include "Mapping.hpp"
#include <boost/config.hpp>
#include <ostream>
#include "mesh/Utils.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

namespace precice::mapping {

Mapping::Mapping(
    Constraint constraint,
    int        dimensions,
    bool       requiresGradientData,
    Type       mappingType)
    : _requiresGradientData(requiresGradientData),
      _constraint(constraint),
      _inputRequirement(MeshRequirement::UNDEFINED),
      _outputRequirement(MeshRequirement::UNDEFINED),
      _input(),
      _output(),
      _dimensions(dimensions),
      _mappingType(mappingType)
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

bool Mapping::isIterative() const
{
  return _mappingType == Type::Iterative;
}

bool Mapping::hasInitialGuess() const
{
  PRECICE_ASSERT(isIterative(), "This mapping isn't iterative, so it cannot have an initial guess.");
  PRECICE_ASSERT(_initialGuess != nullptr, "The last solution wasn't provided.");
  return _initialGuess->size() > 0;
}

const Eigen::VectorXd &Mapping::initialGuess() const
{
  PRECICE_ASSERT(isIterative(), "This mapping isn't iterative, so it cannot compute an initial guess.");
  PRECICE_ASSERT(_initialGuess != nullptr, "The last solution wasn't provided.");
  return *_initialGuess;
}

Eigen::VectorXd &Mapping::initialGuess()
{
  PRECICE_ASSERT(isIterative(), "This mapping isn't iterative.");
  PRECICE_ASSERT(_initialGuess != nullptr, "The last solution wasn't provided.");
  return *_initialGuess;
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

bool Mapping::requiresGradientData() const
{
  return _requiresGradientData;
}

void Mapping::map(int inputDataID, int outputDataID, Eigen::VectorXd &initialGuess)
{
  PRECICE_ASSERT(_initialGuess == nullptr);
  _initialGuess = &initialGuess;
  map(inputDataID, outputDataID);
  _initialGuess = nullptr;
}

void Mapping::map(int inputDataID,
                  int outputDataID)
{
  PRECICE_ASSERT(_hasComputedMapping);
  PRECICE_ASSERT(!isIterative() || _initialGuess != nullptr, "Call the map version with lastSolution");

  PRECICE_ASSERT(input()->getDimensions() == output()->getDimensions(),
                 input()->getDimensions(), output()->getDimensions());
  PRECICE_ASSERT(getDimensions() == output()->getDimensions(),
                 getDimensions(), output()->getDimensions());
  PRECICE_ASSERT(input()->data(inputDataID)->getDimensions() == output()->data(outputDataID)->getDimensions(),
                 input()->data(inputDataID)->getDimensions(), output()->data(outputDataID)->getDimensions());
  PRECICE_ASSERT(input()->data(inputDataID)->values().size() / input()->data(inputDataID)->getDimensions() == static_cast<int>(input()->vertices().size()),
                 input()->data(inputDataID)->values().size(), input()->data(inputDataID)->getDimensions(), input()->vertices().size());
  PRECICE_ASSERT(output()->data(outputDataID)->values().size() / output()->data(outputDataID)->getDimensions() == static_cast<int>(output()->vertices().size()),
                 output()->data(outputDataID)->values().size(), output()->data(outputDataID)->getDimensions(), output()->vertices().size());

  time::Sample sample{input()->data(inputDataID)->getDimensions(),
                      input()->data(inputDataID)->values(),
                      input()->data(inputDataID)->gradients()};
  map(sample, output()->data(outputDataID)->values());
}

void Mapping::map(const time::Sample &input, Eigen::VectorXd &output, Eigen::VectorXd &lastSolution)
{
  PRECICE_ASSERT(_initialGuess == nullptr);
  _initialGuess = &lastSolution;
  map(input, output);
  _initialGuess = nullptr;
}

void Mapping::map(const time::Sample &input, Eigen::VectorXd &output)
{
  PRECICE_ASSERT(_hasComputedMapping);
  PRECICE_ASSERT(!isIterative() || _initialGuess != nullptr, "Call the map version with lastSolution");

  if (hasConstraint(CONSERVATIVE)) {
    mapConservative(input, output);
  } else if (hasConstraint(CONSISTENT)) {
    mapConsistent(input, output);
  } else if (isScaledConsistent()) {
    mapConsistent(input, output);
    scaleConsistentMapping(input.values, output, getConstraint());
  } else {
    PRECICE_UNREACHABLE("Unknown mapping constraint.")
  }
}

void Mapping::scaleConsistentMapping(const Eigen::VectorXd &input, Eigen::VectorXd &output, Mapping::Constraint constraint) const
{
  PRECICE_ASSERT(isScaledConsistent());

  if (input.size() == 0 || output.size() == 0) {
    return;
  }
  PRECICE_ASSERT(input.size() > 0 && output.size() > 0);

  bool            volumeMode = hasConstraint(SCALED_CONSISTENT_VOLUME);
  logging::Logger _log{"mapping::Mapping"};
  // Only serial participant is supported for scale-consistent mapping
  PRECICE_ASSERT((not utils::IntraComm::isPrimary()) and (not utils::IntraComm::isSecondary()));

  // If rank is not empty and do not contain connectivity information, raise error
  int  spaceDimension    = this->input()->getDimensions();
  bool requiresEdges     = (spaceDimension == 2 and !volumeMode);
  bool requiresTriangles = (spaceDimension == 2 and volumeMode) or (spaceDimension == 3 and !volumeMode);
  bool requiresTetra     = (spaceDimension == 3 and volumeMode);

  for (mesh::PtrMesh mesh : {this->input(), this->output()}) {
    if (not mesh->vertices().empty()) {

      PRECICE_CHECK(!(requiresEdges && mesh->edges().empty()), "Edges connectivity information is missing for the mesh \"{}\". "
                                                               "Scaled consistent mapping requires connectivity information.",
                    mesh->getName());

      PRECICE_CHECK(!(requiresTriangles && mesh->triangles().empty()), "Triangles connectivity information is missing for the mesh \"{}\". "
                                                                       "Scaled consistent mapping requires connectivity information.",
                    mesh->getName());

      PRECICE_CHECK(!(requiresTetra && mesh->tetrahedra().empty()), "Tetrahedra connectivity information is missing for the mesh \"{}\". "
                                                                    "Scaled consistent mapping requires connectivity information.",
                    mesh->getName());
    }
  }

  const int valueDimensions = input.size() / this->input()->vertices().size();

  Eigen::VectorXd integralInput;
  Eigen::VectorXd integralOutput;

  // Integral is calculated on each direction separately
  if (!volumeMode) {
    integralInput  = mesh::integrateSurface(this->input(), input);
    integralOutput = mesh::integrateSurface(this->output(), output);
  } else {
    integralInput  = mesh::integrateVolume(this->input(), input);
    integralOutput = mesh::integrateVolume(this->output(), output);
  }

  // Create reshape the output values vector to matrix
  Eigen::Map<Eigen::MatrixXd> outputValuesMatrix(output.data(), valueDimensions, output.size() / valueDimensions);

  // Scale in each direction
  Eigen::VectorXd scalingFactor = integralInput.array() / integralOutput.array();
  PRECICE_DEBUG("Scaling factor in scaled-consistent mapping: {}", scalingFactor);
  outputValuesMatrix.array().colwise() *= scalingFactor.array();
} // namespace mapping

bool Mapping::hasConstraint(const Constraint &constraint) const
{
  return (getConstraint() == constraint);
}

bool Mapping::hasComputedMapping() const
{
  return _hasComputedMapping;
}

bool Mapping::isScaledConsistent() const
{
  return (hasConstraint(SCALED_CONSISTENT_SURFACE) || hasConstraint(SCALED_CONSISTENT_VOLUME));
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

} // namespace precice::mapping
