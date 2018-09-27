#include <Eigen/Core>
#include <mesh/RTree.hpp>
#include "NearestProjectionMapping.hpp"
#include "query/FindClosest.hpp"
#include "utils/EventTimings.hpp"

namespace precice
{
namespace mapping
{

NearestProjectionMapping:: NearestProjectionMapping
(
  Constraint constraint,
  int        dimensions)
:
  Mapping(constraint, dimensions)
{
  if (constraint == CONSISTENT) {
    setInputRequirement(FULL);
    setOutputRequirement(VERTEX);
  } else {
    assertion(constraint == CONSERVATIVE, constraint);
    setInputRequirement(VERTEX);
    setOutputRequirement(FULL);
  }
}

void NearestProjectionMapping::computeMapping()
{
  TRACE(input()->vertices().size(), output()->vertices().size());
  precice::utils::Event e("map.np.computeMapping.From" + input()->getName() + "To" + output()->getName());

  if (getConstraint() == CONSISTENT){
    DEBUG("Compute consistent mapping");
    auto        rtree     = indexMesh(*input());
    const auto &oVertices = output()->vertices();
    _weights.resize(oVertices.size());
    for (size_t i = 0; i < oVertices.size(); i++) {
      const Eigen::VectorXd &coords = oVertices[i].getCoords();
      // Search for the output vertex inside the input mesh
      rtree.query(boost::geometry::index::nearest(coords, 1),
                  boost::make_function_output_iterator([&](const mesh::PrimitiveRTree::value_type &pnearest) {
                      using query::generateInterpolationElements;
                      using mesh::Primitive;
                    const auto& nearest = pnearest.second;
                    // fill the weights
                    switch (nearest.type) {
                    case (Primitive::Vertex):
                      _weights[i] = generateInterpolationElements(oVertices[i],
                                                                  input()->vertices()[nearest.index]);
                      break;
                    case (Primitive::Edge):
                      _weights[i] = generateInterpolationElements(oVertices[i],
                                                                  input()->edges()[nearest.index]);
                      break;
                    case (Primitive::Triangle):
                      _weights[i] = generateInterpolationElements(oVertices[i],
                                                                  input()->triangles()[nearest.index]);
                      break;
                    case (Primitive::Quad):
                      _weights[i] = generateInterpolationElements(oVertices[i],
                                                                  input()->quads()[nearest.index]);
                      break;
                    }
                    CHECK(!_weights[i].empty(),
                          "No interpolation elements for current vertex!");
                  }));
    }
    assertion(std::none_of(_weights.cbegin(), _weights.cend(), [](const query::InterpolationElements &elements) {
              return elements.empty();
            }),
            "The mapping is incomplete as there are vertices with no interpolation elements assigned to them.");
  } else {
    assertion(getConstraint() == CONSERVATIVE, getConstraint());
    DEBUG("Compute conservative mapping");
    auto        rtree     = indexMesh(*output());
    const auto &iVertices = input()->vertices();
    _weights.resize(iVertices.size());
    for (size_t i = 0; i < iVertices.size(); i++) {
      const Eigen::VectorXd &coords = iVertices[i].getCoords();
      // Search for the output vertex inside the input mesh
      rtree.query(boost::geometry::index::nearest(coords, 1),
                  boost::make_function_output_iterator([&](const mesh::PrimitiveRTree::value_type &pnearest) {
                      using query::generateInterpolationElements;
                      using mesh::Primitive;
                    const auto& nearest = pnearest.second;
                    switch (nearest.type) {
                    case (Primitive::Vertex):
                      _weights[i] = generateInterpolationElements(iVertices[i],
                                                                  output()->vertices()[nearest.index]);
                      break;
                      case (Primitive::Edge):
                      _weights[i] = generateInterpolationElements(iVertices[i],
                                                                  output()->edges()[nearest.index]);
                      break;
                      case (Primitive::Triangle):
                      _weights[i] = generateInterpolationElements(iVertices[i],
                                                                  output()->triangles()[nearest.index]);
                      break;
                      case (Primitive::Quad):
                      _weights[i] = generateInterpolationElements(iVertices[i],
                                                                  output()->quads()[nearest.index]);
                      break;
                    }
                    CHECK(!_weights[i].empty(),
                          "No interpolation elements for current vertex!");
                  }));
    }
    assertion(std::none_of(_weights.cbegin(), _weights.cend(), [](const query::InterpolationElements &elements) {
              return elements.empty();
            }),
            "The mapping is incomplete as there are vertices with no interpolation elements assigned to them.");
  }
  _hasComputedMapping = true;
}

bool NearestProjectionMapping::hasComputedMapping() const
{
  return _hasComputedMapping;
}

void NearestProjectionMapping::clear()
{
  TRACE();
  _weights.clear();
  _hasComputedMapping = false;
}

void NearestProjectionMapping::map(
    int inputDataID,
    int outputDataID)
{
  TRACE(inputDataID, outputDataID);
  precice::utils::Event e("map.np.mapData.From" + input()->getName() + "To" + output()->getName());

  mesh::PtrData inData = input()->data(inputDataID);
  mesh::PtrData outData = output()->data(outputDataID);
  const Eigen::VectorXd& inValues = inData->values();
  Eigen::VectorXd& outValues = outData->values();
  int dimensions = inData->getDimensions();
  assertion(dimensions == outData->getDimensions());

  if (getConstraint() == CONSISTENT) {
    DEBUG("Map consistent");
    assertion(_weights.size() == output()->vertices().size(),
              _weights.size(), output()->vertices().size());
    for (size_t i = 0; i < output()->vertices().size(); i++) {
      InterpolationElements &elems     = _weights[i];
      size_t                 outOffset = i * dimensions;
      for (query::InterpolationElement &elem : elems) {
        size_t inOffset = (size_t) elem.element->getID() * dimensions;
        for (int dim = 0; dim < dimensions; dim++) {
          assertion(outOffset + dim < (size_t) outValues.size());
          assertion(inOffset + dim < (size_t) inValues.size());
          outValues(outOffset + dim) += elem.weight * inValues(inOffset + dim);
        }
      }
    }
  } else {
    assertion(getConstraint() == CONSERVATIVE, getConstraint());
    DEBUG("Map conservative");
    assertion(_weights.size() == input()->vertices().size(),
              _weights.size(), input()->vertices().size());
    for (size_t i = 0; i < input()->vertices().size(); i++) {
      size_t                 inOffset = i * dimensions;
      InterpolationElements &elems    = _weights[i];
      for (query::InterpolationElement &elem : elems) {
        size_t outOffset = (size_t) elem.element->getID() * dimensions;
        for (int dim = 0; dim < dimensions; dim++) {
          assertion(outOffset + dim < (size_t) outValues.size());
          assertion(inOffset + dim < (size_t) inValues.size());
          outValues(outOffset + dim) += elem.weight * inValues(inOffset + dim);
        }
      }
    }
  }
}

void NearestProjectionMapping::tagMeshFirstRound()
{
  TRACE();

  computeMapping();

  if (getConstraint() == CONSISTENT) {
    for (mesh::Vertex &v : input()->vertices()) {
      for (size_t i = 0; i < output()->vertices().size(); i++) {
        const InterpolationElements &elems = _weights[i];
        for (const query::InterpolationElement &elem : elems) {
          if (elem.element->getID() == v.getID() && elem.weight != 0.0) {
            v.tag();
          }
        }
      }
    }
  } else {
    assertion(getConstraint() == CONSERVATIVE, getConstraint());
    for (mesh::Vertex &v : output()->vertices()) {
      for (size_t i = 0; i < input()->vertices().size(); i++) {
        const InterpolationElements &elems = _weights[i];
        for (const query::InterpolationElement &elem : elems) {
          if (elem.element->getID() == v.getID() && elem.weight != 0.0) {
            v.tag();
          }
        }
      }
    }
  }

  clear();
}

void NearestProjectionMapping::tagMeshSecondRound()
{
  TRACE();
  // for NP mapping no operation needed here
}

} // namespace mapping
} // namespace precice
