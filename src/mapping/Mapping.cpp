#include "Mapping.hpp"
#include <boost/config.hpp>
#include <ostream>
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
  std::vector<double> integralInput(valueDimensions);
  std::vector<double> integralOutput(valueDimensions);

  // Initialize integral values
  std::fill(integralInput.begin(), integralInput.end(), 0.0);
  std::fill(integralOutput.begin(), integralOutput.end(), 0.0);

  // Compute integral on input and output meshes
  if (meshDimensions == 2) {
    // Calculate integral on the input mesh
    for (const auto &edge : input()->edges()) {
      double area    = edge.getLength();
      int    vertex1 = edge.vertex(0).getID() * valueDimensions;
      int    vertex2 = edge.vertex(1).getID() * valueDimensions;
      // Input mesh may have overlaps, we only need to include the values on the owned vertices
      if (edge.vertex(0).isOwner() and edge.vertex(1).isOwner()) {
        for (int dim = 0; dim < valueDimensions; ++dim) {
          integralInput.at(dim) += 0.5 * area * (inputValues(vertex1 + dim) + inputValues(vertex2 + dim));
        }
      }
    }
    // Calculate on the output mesh
    for (const auto &edge : output()->edges()) {
      double area    = edge.getLength();
      int    vertex1 = edge.vertex(0).getID() * valueDimensions;
      int    vertex2 = edge.vertex(1).getID() * valueDimensions;
      for (int dim = 0; dim < valueDimensions; ++dim) {
        integralOutput.at(dim) += 0.5 * area * (outputValues(vertex1 + dim) + outputValues(vertex2 + dim));
      }
    }
  } else { // 3D
    // For 3D case, only edge connectivity is not enough. Also the face connectivity is necessary.
    PRECICE_ASSERT(not input()->triangles().empty() or not output()->triangles().empty());
    // Calculate integral on the input mesh
    for (const auto &face : input()->triangles()) {
      double area    = face.getArea();
      int    vertex1 = face.vertex(0).getID() * valueDimensions;
      int    vertex2 = face.vertex(1).getID() * valueDimensions;
      int    vertex3 = face.vertex(2).getID() * valueDimensions;
      // Input mesh may have overlaps, we only need to include the values on the owned vertices
      if (face.vertex(0).isOwner() and face.vertex(1).isOwner() and face.vertex(2).isOwner()) {
        for (int dim = 0; dim < valueDimensions; ++dim) {
          integralInput.at(dim) += (area / 3.0) * (inputValues(vertex1 + dim) + inputValues(vertex2 + dim) + inputValues(vertex3 + dim));
        }
      }
    }
    // Calculate integral on the output mesh
    for (const auto &face : output()->triangles()) {
      double area    = face.getArea();
      int    vertex1 = face.vertex(0).getID() * valueDimensions;
      int    vertex2 = face.vertex(1).getID() * valueDimensions;
      int    vertex3 = face.vertex(2).getID() * valueDimensions;
      for (int dim = 0; dim < valueDimensions; ++dim) {
        integralOutput.at(dim) += (area / 3.0) * (outputValues(vertex1 + dim) + outputValues(vertex2 + dim) + outputValues(vertex3 + dim));
      }
    }
  }

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
