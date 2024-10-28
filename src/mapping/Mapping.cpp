#include "Mapping.hpp"
#include <boost/config.hpp>
#include <ostream>
#include "math/differences.hpp"
#include "mesh/Utils.hpp"
#include "profiling/Event.hpp"
#include "utils/IntraComm.hpp"
#include "utils/assertion.hpp"

namespace precice::mapping {

Mapping::Mapping(
    Constraint              constraint,
    int                     dimensions,
    bool                    requiresGradientData,
    InitialGuessRequirement mappingType)
    : _requiresGradientData(requiresGradientData),
      _constraint(constraint),
      _inputRequirement(MeshRequirement::UNDEFINED),
      _outputRequirement(MeshRequirement::UNDEFINED),
      _input(),
      _output(),
      _dimensions(dimensions),
      _initialGuessRequirement(mappingType)
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

bool Mapping::requiresInitialGuess() const
{
  return _initialGuessRequirement == InitialGuessRequirement::Required;
}

bool Mapping::hasInitialGuess() const
{
  PRECICE_ASSERT(requiresInitialGuess(), "This mapping isn't iterative, so it cannot have an initial guess.");
  PRECICE_ASSERT(_initialGuess != nullptr, "The last solution wasn't provided.");
  return _initialGuess->size() > 0;
}

const Eigen::VectorXd &Mapping::initialGuess() const
{
  PRECICE_ASSERT(requiresInitialGuess(), "This mapping isn't iterative, so it doesn't have an initial guess.");
  PRECICE_ASSERT(_initialGuess != nullptr, "The last solution wasn't provided.");
  return *_initialGuess;
}

Eigen::VectorXd &Mapping::initialGuess()
{
  PRECICE_ASSERT(requiresInitialGuess(), "This mapping isn't iterative, so it doesn't have an initial guess.");
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
  PRECICE_ASSERT(!requiresInitialGuess() || _initialGuess != nullptr, "Call the map version with lastSolution");

  PRECICE_ASSERT(input()->getDimensions() == output()->getDimensions(),
                 input()->getDimensions(), output()->getDimensions());
  PRECICE_ASSERT(getDimensions() == output()->getDimensions(),
                 getDimensions(), output()->getDimensions());
  PRECICE_ASSERT(input()->data(inputDataID)->getDimensions() == output()->data(outputDataID)->getDimensions(),
                 input()->data(inputDataID)->getDimensions(), output()->data(outputDataID)->getDimensions());
  PRECICE_ASSERT(input()->data(inputDataID)->values().size() / input()->data(inputDataID)->getDimensions() == static_cast<int>(input()->nVertices()),
                 input()->data(inputDataID)->values().size(), input()->data(inputDataID)->getDimensions(), input()->nVertices());
  PRECICE_ASSERT(output()->data(outputDataID)->values().size() / output()->data(outputDataID)->getDimensions() == static_cast<int>(output()->nVertices()),
                 output()->data(outputDataID)->values().size(), output()->data(outputDataID)->getDimensions(), output()->nVertices());

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
  PRECICE_ASSERT(!requiresInitialGuess() || _initialGuess != nullptr, "Call the map version with lastSolution");

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
    if (not mesh->empty()) {

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

  const int valueDimensions = input.size() / this->input()->nVertices();

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
  // We cannot handle the case with zero output data and non-zero input data.
  // To fulfill the constraint, we would need to scale the output data in such a way, that the integral sum of the input is preserved.
  // That's not possible using a constant scaling factor as the output sum will always be zero. Here, we return 1 and emit a warning afterwards.
  const Eigen::VectorXd scalingFactor = integralInput.binaryExpr(integralOutput, [](double lhs, double rhs) { return (rhs == 0.0) ? 1.0 : (lhs / rhs); });
  PRECICE_DEBUG("Scaling factor in scaled-consistent mapping: {}", scalingFactor);

  PRECICE_DEBUG("Scaling factor in scaled-consistent mapping: {}", scalingFactor);
  outputValuesMatrix.array().colwise() *= scalingFactor.array();

  // check whether the constraint is fulfilled
  for (Eigen::Index i = 0; i < scalingFactor.size(); ++i) {
    double consistency = scalingFactor[i] * integralOutput[i] - integralInput[i];
    PRECICE_WARN_IF(
        math::greater(std::abs(consistency), 0.0),
        "Failed to fulfill consistency constraint of component {} for scaled-consistent mapping from mesh \"{}\" to mesh \"{}\". Consistency difference between input and scaled output is \"{}\".", i, this->input()->getName(), this->output()->getName(), consistency);
  }
}

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

bool Mapping::isIndirectMapping() const
{
  return _isIndirect;
}

void Mapping::updateMappingDataCache(MappingDataCache &cache, Eigen::VectorXd &in)
{
  precice::profiling::Event e("map.updateCache.From" + input()->getName());
  cache.inData = in;
}

void Mapping::writeConservativeAt(::precice::span<const double> coordinates, Eigen::Map<const Eigen::MatrixXd> &source, Eigen::Map<Eigen::MatrixXd> &target)
{
  PRECICE_ASSERT(false, "Not implemented");
}

void Mapping::evaluateMappingDataCacheAt(::precice::span<const double> coordinates, const MappingDataCache &cache, ::precice::span<double> values)
{
  PRECICE_ASSERT(false, "Not implemented");
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
