#include "CoarseGrainingMapping.hpp"

#include <boost/container/flat_set.hpp>
#include "logging/LogMacros.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/RadialBasisFctSolver.hpp"
#include "math/math.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/Vertex.hpp"
#include "profiling/Event.hpp"
#include "utils/IntraComm.hpp"
#include "utils/Parallel.hpp"

namespace precice::mapping {

namespace impl {
class LucyKernelFunction {
public:
  LucyKernelFunction(short grainDim, double functionRadius)
      : _grainDim(grainDim), _c(functionRadius) {}

  double getFunctionRadius() const
  {
    return _c;
  }

  double evaluate(double r) const
  {
    if (r > _c) {
      return 0;
    }

    constexpr double pi     = 3.1415926535897931;
    const double     ratio  = r / _c;
    double           res    = 0;
    double           factor = 0;

    // TODO: We could directly store inv c in the class itself
    switch (_grainDim) {
    case 1:
      factor = 5. * (1. / (3 * _c));
      res    = factor * ((1 + ratio) * math::pow_int<3>(1 - ratio));
      break;
    case 2:
      factor = 10. * (1. / (pi * math::pow_int<2>(_c)));
      res    = factor * (math::pow_int<3>(1 - ratio));
      break;
    case 3:
      factor = 105. * (1. / (16 * pi * math::pow_int<3>(_c)));
      res    = factor * ((1 + 3 * ratio) * math::pow_int<3>(1 - ratio));
      break;
    default:
      PRECICE_UNREACHABLE("Unknown dimension");
      break;
    }
    return res;
  }

private:
  const short  _grainDim{};
  const double _c{};
};

} // namespace impl

CoarseGrainingMapping::CoarseGrainingMapping(
    Constraint constraint,
    int        meshDim,
    int        grainDim,
    double     functionRadius)
    : Mapping(constraint, meshDim, /*requiresGradientData*/ false, Mapping::InitialGuessRequirement::None)
{
  PRECICE_CHECK(functionRadius > 0, "Function radius must be greater zero.");
  _lucyFunction = std::make_unique<impl::LucyKernelFunction>(static_cast<short>(grainDim), functionRadius);
}

void CoarseGrainingMapping::mapConsistentAt(const Eigen::Ref<const Eigen::MatrixXd> &coordinates, const impl::MappingDataCache &cache, Eigen::Ref<Eigen::MatrixXd> values)
{
  PRECICE_CHECK(false, "consistent constraint is not implemented.");
}

void CoarseGrainingMapping::mapConservativeAt(const Eigen::Ref<const Eigen::MatrixXd> &coordinates, const Eigen::Ref<const Eigen::MatrixXd> &source, impl::MappingDataCache &cache, Eigen::Ref<Eigen::MatrixXd> target)
{
  precice::profiling::Event e("map.cg.mapConservativeAt.From" + input()->getName());
  auto                     &index = output()->index();
  for (Eigen::Index i = 0; i < coordinates.cols(); ++i) {
    mesh::Vertex src{coordinates.col(i), -1};
    auto         dest = index.getVerticesInsideBox(src, _lucyFunction->getFunctionRadius());

    for (const auto &d : dest) {
      const auto &dst   = output()->vertex(d).rawCoords();
      auto        dist  = computeSquaredDifference(dst, src.rawCoords());
      auto        coeff = _lucyFunction->evaluate(dist);
      target.col(d) += coeff * source.col(i);
    }
  }
}

void CoarseGrainingMapping::computeMapping()
{
  PRECICE_TRACE(input()->nVertices());

  PRECICE_ASSERT(input().get() != nullptr);
  PRECICE_ASSERT(output().get() != nullptr);

  const std::string         baseEvent = "map.cg.computeMapping.From" + input()->getName() + "To" + output()->getName();
  precice::profiling::Event e(baseEvent, profiling::Synchronize);

  // Setup Direction of Mapping
  PRECICE_CHECK(hasConstraint(CONSERVATIVE), "Only conservative constraints are implemented.");
  mesh::PtrMesh in  = input();
  mesh::PtrMesh out = output();

  // Set up of output arrays
  const size_t verticesSize   = in->nVertices();
  const auto  &sourceVertices = in->vertices();
  // _vertexIndices.resize(verticesSize);

  // Needed for error calculations
  // utils::statistics::DistanceAccumulator distanceStatistics;

  // auto &index = out->index();
  // for (size_t i = 0; i < verticesSize; ++i) {
  //   const auto &sourceCoords  = sourceVertices[i].getCoords();
  //   const auto  matchedVertex = index.getClosestVertex(sourceCoords);
  //   _vertexIndices[i]         = matchedVertex.index;

  //   // Compute distance between input and output vertiex for the stats
  //   const auto &matchCoords = out->vertex(matchedVertex.index).getCoords();
  //   auto        distance    = (sourceCoords - matchCoords).norm();
  //   distanceStatistics(distance);
  // }
  _hasComputedMapping = true;
}

void CoarseGrainingMapping::mapConservative(const time::Sample &inData, Eigen::VectorXd &outData)
{
  PRECICE_CHECK(false, "only just-in-time-mapping variant is implemented.");
}

void CoarseGrainingMapping::mapConsistent(const time::Sample &inData, Eigen::VectorXd &outData)
{
  PRECICE_CHECK(false, "only just-in-time-mapping conservative variant is implemented.");
}

void CoarseGrainingMapping::clear()
{
  PRECICE_TRACE();
  // _vertexIndices.clear();
  _hasComputedMapping = false;

  if (getConstraint() == CONSISTENT) {
    input()->index().clear();
  } else {
    output()->index().clear();
  }
}

std::string CoarseGrainingMapping::getName() const
{
  return "coarse-graining";
}

void CoarseGrainingMapping::tagMeshFirstRound()
{
  PRECICE_TRACE();
  precice::profiling::Event e("map.cg.tagMeshFirstRound.From" + input()->getName() + "To" + output()->getName(), profiling::Synchronize);

  // parallel partitioning for just-in-time mapping:
  if (this->isJustInTimeMapping()) {
    // in the usual case, we make use of the indexSet, which is pre-computed from the mapping
    // for the just-in-time mapping, we can't do that since we don't have the output (local) mesh
    // what we would need to do in theory for a perfect partitioning:
    // find all nearest-neighbors at the 'boundary' of the access region, which would require an
    // infinite fine sampling of output mesh nodes to be used in the computeMapping below
    // for now, we simply tag everything and move on. The remote mesh is here already filtered
    // through the geometric filter setting.
    //
    // Depending on the mapping constraint, one of these tagAll calls will do nothing, as the vertex
    // set of the mesh is empty. From a practical point of view, we only need to apply the
    // tagging to one of the meshes (the remote one). But calling it on both sides reliefs us from any
    // conditional code.
    output()->tagAll();
    input()->tagAll();
    return;
  }

  // computeMapping();

  // // Lookup table of all indices used in the mapping
  // const boost::container::flat_set<int> indexSet(_vertexIndices.begin(), _vertexIndices.end());

  // // Get the source mesh depending on the constraint
  // const mesh::PtrMesh &source = hasConstraint(CONSERVATIVE) ? output() : input();

  // // Tag all vertices used in the mapping
  // for (mesh::Vertex &v : source->vertices()) {
  //   if (indexSet.count(v.getID()) != 0) {
  //     v.tag();
  //   }
  // }

  // clear();
}

void CoarseGrainingMapping::tagMeshSecondRound()
{
  PRECICE_TRACE();
  // for CG mapping no operation needed here
}

} // namespace precice::mapping
