#pragma once

#include <mapping/MathHelper.hpp>
#include <mapping/impl/BasisFunctions.hpp>
#include <mesh/Mesh.hpp>

namespace precice::mapping {

struct Sample {
  double pos;
  double error;

  Sample(double pos, double error)
      : pos(pos), error(error)
  {
  }

  Sample()
      : pos(std::numeric_limits<double>::quiet_NaN()), error(std::numeric_limits<double>::quiet_NaN())
  {
  }
};

template <typename Solver>
class RBFParameterTuner {

  mutable logging::Logger _log{"mapping::RBFParameterTuner"};

  using RBF_T = typename Solver::BASIS_FUNCTION_T; // TODO: better?

protected:
  Eigen::Index _inSize;

  double _lowerBound;
  double _upperBound;
  bool   _lastSampleWasOptimum;

  static double estimateMeshResolution(const mesh::Mesh &inputMesh);

public:
  virtual ~RBFParameterTuner() = default;
  explicit RBFParameterTuner(const mesh::Mesh &inputMesh);

  virtual std::tuple<double, double> optimize(const Solver &solver, const Eigen::VectorXd &inputData);

  bool lastSampleWasOptimum() const;
};

template <typename Solver>
RBFParameterTuner<Solver>::RBFParameterTuner(const mesh::Mesh &inputMesh)
{
  constexpr bool radiusRBF = RadiusInitialization<typename Solver::BASIS_FUNCTION_T>::isAvailable();
  PRECICE_ASSERT(radiusRBF, "RBF is not supported by this optimizer, as it does not accept a support-radius."
                            "Currently supported: Compactly supported RBFs and Gaussians.");

  _lowerBound = estimateMeshResolution(inputMesh);
  _inSize     = inputMesh.nVertices();
  _upperBound = std::numeric_limits<double>::quiet_NaN();

  _lastSampleWasOptimum = false;
}

template <typename Solver>
std::tuple<double, double> RBFParameterTuner<Solver>::optimize(const Solver &solver, const Eigen::VectorXd &inputData)
{
  PRECICE_ASSERT(false, "Not implemented.");
  constexpr double nan = std::numeric_limits<double>::quiet_NaN();
  return {nan, nan};
}

template <typename Solver>
double RBFParameterTuner<Solver>::estimateMeshResolution(const mesh::Mesh &inputMesh)
{
  constexpr int sampleSize = 3;

  const size_t       i0 = inputMesh.vertices().size() / 2;
  const mesh::Vertex x0 = inputMesh.vertices().at(i0);

  const std::vector<int> matches = inputMesh.index().getClosestVertices(x0.getCoords(), sampleSize);

  double h = 0;
  for (int i = 0; i < sampleSize; i++) {
    const mesh::Vertex xi = inputMesh.vertices().at(matches.at(i));
    h += std::sqrt(utils::computeSquaredDifference(xi.rawCoords(), x0.rawCoords()));
  }
  return h / sampleSize;
}

template <typename Solver>
bool RBFParameterTuner<Solver>::lastSampleWasOptimum() const
{
  return _lastSampleWasOptimum;
}

} // namespace precice::mapping
