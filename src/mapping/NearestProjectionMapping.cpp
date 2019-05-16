#include "NearestProjectionMapping.hpp"
#include "query/FindClosest.hpp"
#include <Eigen/Core>
#include "utils/Event.hpp"
#include "mesh/RTree.hpp"
#include <stdexcept>

namespace precice {
extern bool syncMode;

namespace mapping {

NearestProjectionMapping:: NearestProjectionMapping
(
  Constraint constraint,
  int        dimensions)
:
  Mapping(constraint, dimensions)
{
  if (constraint == CONSISTENT){
    setInputRequirement(Mapping::MeshRequirement::FULL);
    setOutputRequirement(Mapping::MeshRequirement::VERTEX);
  }
  else {
    assertion(constraint == CONSERVATIVE, constraint);
    setInputRequirement(Mapping::MeshRequirement::VERTEX);
    setOutputRequirement(Mapping::MeshRequirement::FULL);
  }
}

namespace {
class InterpolationElementsGenerator {
public:
  InterpolationElementsGenerator(const mesh::Mesh &mesh)
      : _mesh(mesh){};

  query::InterpolationElements operator()(const mesh::Vertex& pos, mesh::PrimitiveIndex pi) const
  {
    using query::generateInterpolationElements;
    using mesh::Primitive;
    const auto idx = pi.index;
    switch (pi.type) {
    case (Primitive::Vertex):
      return generateInterpolationElements(pos, _mesh.vertices()[idx]);
    case (Primitive::Edge):
      return generateInterpolationElements(pos, _mesh.edges()[idx]);
    case (Primitive::Triangle):
      return generateInterpolationElements(pos, _mesh.triangles()[idx]);
    case (Primitive::Quad):
      return generateInterpolationElements(pos, _mesh.quads()[idx]);
    default:
      throw std::invalid_argument{"Primitve is unknown"};
    }
  }

private:
  const mesh::Mesh &_mesh;
};
} // namespace

void NearestProjectionMapping:: computeMapping()
{
  TRACE()
  precice::utils::Event e("map.np.computeMapping.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);
  if (getConstraint() == CONSISTENT){
      computeMappingConsistent();
  } else {
      computeMappingConservative();
  }
}

void NearestProjectionMapping:: computeMappingConsistent()
{
    TRACE(input()->vertices().size(), output()->vertices().size());
    DEBUG("Compute consistent mapping");
    assertion(getConstraint() == CONSISTENT, getConstraint());

    if (getDimensions() == 2) {
        auto rtree = mesh::rtree::getEdgeRTree(input());
        const auto &oVertices = output()->vertices();
        const auto &iEdges    = input()->edges();
        CHECK(oVertices.empty() || !iEdges.empty(), "Mesh \"" << input()->getName() << "\" does not contain Edges to project onto.");
        _weights.resize(oVertices.size());
        for (size_t i = 0; i < oVertices.size(); i++) {
            const Eigen::VectorXd &coords = oVertices[i].getCoords();
            // Search for the output vertex inside the input mesh
            rtree->query(boost::geometry::index::nearest(coords, 1),
                     boost::make_function_output_iterator([&](size_t eid) {
                         auto& weights = _weights[i];
                         weights = query::generateInterpolationElements(oVertices[i], iEdges[eid]);
                         CHECK(!weights.empty(), "No interpolation elements for current vertex!");
                         if (std::any_of(weights.begin(), weights.end(), [](const query::InterpolationElement& elem){return elem.weight < 0;})) {
                            WARN("Mapping \"" << input()->getName() << "\" contains vertex (" << oVertices[i] << ") which has negative weights indicating non-matching meshes!");
                         }
                    }));
        }
    } else {
        auto rtree = mesh::rtree::getTriangleRTree(input());
        const auto &oVertices  = output()->vertices();
        const auto &iTriangles = input()->triangles();
        CHECK(oVertices.empty() || !iTriangles.empty(), "Mesh \"" << input()->getName() << "\" does not contain Triangles to project onto.");
        _weights.resize(oVertices.size());
        using IndexType = typename mesh::rtree::triangle_traits::IndexType;
        std::vector<IndexType> matches;
        matches.reserve(10);
        for (size_t i = 0; i < oVertices.size(); i++) {
            const Eigen::VectorXd &coords = oVertices[i].getCoords();
            // Search for the output vertex inside the input mesh
            rtree->query(boost::geometry::index::nearest(coords, 10), std::back_inserter(matches));
            auto closestId = std::min_element(matches.begin(), matches.end(), 
                    [&](const IndexType& lhs, const IndexType& rhs) {
                    using boost::geometry::comparable_distance;
                    auto distLhs = comparable_distance(
                            coords, iTriangles[lhs.second]
                            );
                    auto distRhs = comparable_distance(
                            coords, iTriangles[lhs.second]
                            );
                    return distLhs < distRhs;
                    })->second;
             matches.clear();
             auto& weights = _weights[i];
             weights = query::generateInterpolationElements(oVertices[i], iTriangles[closestId]);
             CHECK(!weights.empty(), "No interpolation elements for current vertex!");
             if (std::any_of(weights.begin(), weights.end(), [](const query::InterpolationElement& elem){return elem.weight < 0;})) {
                WARN("Mapping \"" << input()->getName() << "\" contains vertex (" << oVertices[i] << ") which has negative weights indicating non-matching meshes!");
             }
        }
    }
    _hasComputedMapping = true;
}

void NearestProjectionMapping:: computeMappingConservative()
{
    TRACE(input()->vertices().size(), output()->vertices().size());
    DEBUG("Compute conservative mapping");
    assertion(getConstraint() == CONSERVATIVE, getConstraint());

    if (getDimensions() == 2) {
        auto rtree = mesh::rtree::getEdgeRTree(output());
        const auto &iVertices = input()->vertices();
        const auto &oEdges    = output()->edges();
        CHECK(iVertices.empty() || !oEdges.empty(), "Mesh \"" << output()->getName() << "\" does not contain Edges to project onto.");
        _weights.resize(iVertices.size());
        for (size_t i = 0; i < iVertices.size(); i++) {
            const Eigen::VectorXd &coords = iVertices[i].getCoords();
            // Search for the output vertex inside the input mesh
            rtree->query(boost::geometry::index::nearest(coords, 1),
                     boost::make_function_output_iterator([&](size_t eid) {
                         auto& weights = _weights[i];
                         weights = query::generateInterpolationElements(iVertices[i], oEdges[eid]);
                         CHECK(!weights.empty(), "No interpolation elements for current vertex!");
                         if (std::any_of(weights.begin(), weights.end(), [](const query::InterpolationElement& elem){return elem.weight < 0;})) {
                            WARN("Mesh \"" << output()->getName() << "\" contains vertex (" << iVertices[i] << ") which has negative weights indicating non-matching meshes!");
                         }
                    }));
        }
    } else {
        auto rtree = mesh::rtree::getTriangleRTree(input());
        const auto &iVertices  = input()->vertices();
        const auto &oTriangles = output()->triangles();
        CHECK(iVertices.empty() || !oTriangles.empty(), "Mesh \"" << output()->getName() << "\" does not contain Triangles to project onto.");
        _weights.resize(iVertices.size());
        using IndexType = typename mesh::rtree::triangle_traits::IndexType;
        std::vector<IndexType> matches;
        for (size_t i = 0; i < iVertices.size(); i++) {
            const Eigen::VectorXd &coords = iVertices[i].getCoords();
            // Search for the output vertex inside the input mesh
            rtree->query(boost::geometry::index::nearest(coords, 10), std::back_inserter(matches));
            auto closestId = std::min_element(matches.begin(), matches.end(), 
                    [&](const IndexType& lhs, const IndexType& rhs) {
                    using boost::geometry::comparable_distance;
                    auto distLhs = comparable_distance(
                            coords, oTriangles[lhs.second]
                            );
                    auto distRhs = comparable_distance(
                            coords, oTriangles[lhs.second]
                            );
                    return distLhs < distRhs;
                    })->second;
             matches.clear();
             auto& weights = _weights[i];
             weights = query::generateInterpolationElements(iVertices[i], oTriangles[closestId]);
             CHECK(!weights.empty(), "No interpolation elements for current vertex!");
             if (std::any_of(weights.begin(), weights.end(), [](const query::InterpolationElement& elem){return elem.weight < 0;})) {
                 WARN("Mesh \"" << output()->getName() << "\" contains vertex (" << iVertices[i] << ") which has negative weights indicating non-matching meshes!");
             }
        }
    }
    _hasComputedMapping = true;
}



#if 0
{
    TRACE(input()->vertices().size(), output()->vertices().size());
    assertion(getConstraint() == CONSERVATIVE, getConstraint());
    auto        rtree     = indexMesh(*output());
    InterpolationElementsGenerator gen(*output());
    const auto &iVertices = input()->vertices();
    _weights.resize(iVertices.size());
    for (size_t i = 0; i < iVertices.size(); i++) {
        const Eigen::VectorXd &coords = iVertices[i].getCoords();
        // Search for the output vertex inside the input mesh
        rtree->query(boost::geometry::index::nearest(coords, 1),
                boost::make_function_output_iterator([&](const mesh::PrimitiveRTree::value_type &pnearest) {
                    using query::generateInterpolationElements;
                    using mesh::Primitive;
                    const auto& nearest = pnearest.second;
                    auto& weights = _weights[i];
                    weights = gen(iVertices[i], nearest);
                    CHECK(!weights.empty(),
                            "No interpolation elements for current vertex!");
                    if(std::any_of(weights.begin(), weights.end(), [](const query::InterpolationElement& elem){return elem.weight < 0;})) {
                    WARN("Mapping \"" << input()->getName() << "\" contains vertex (" << iVertices[i] << ") which has negative weights indicating non-matching meshes!");
                    }
                    }));
    }
    assertion(std::none_of(_weights.cbegin(), _weights.cend(), [](const query::InterpolationElements &elements) {
                return elements.empty();
                }),
            "The mapping is incomplete as there are vertices with no interpolation elements assigned to them.");
    _hasComputedMapping = true;
}
#endif

bool NearestProjectionMapping:: hasComputedMapping() const
{
  return _hasComputedMapping;
}

void NearestProjectionMapping:: clear()
{
  TRACE();
  _weights.clear();
  _hasComputedMapping = false;
  if (getConstraint() == CONSISTENT){
    mesh::rtree::clear(*input()); 
  } else {
    mesh::rtree::clear(*output()); 
  }
}

void NearestProjectionMapping:: map
(
  int inputDataID,
  int outputDataID )
{
  TRACE(inputDataID, outputDataID);

  precice::utils::Event e("map.np.mapData.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);

  mesh::PtrData inData = input()->data(inputDataID);
  mesh::PtrData outData = output()->data(outputDataID);
  const Eigen::VectorXd& inValues = inData->values();
  Eigen::VectorXd& outValues = outData->values();

  int dimensions = inData->getDimensions();
  assertion(dimensions == outData->getDimensions());

  if (getConstraint() == CONSISTENT){
    DEBUG("Map consistent");
    assertion(_weights.size() == output()->vertices().size(),
               _weights.size(), output()->vertices().size());
    for (size_t i=0; i < output()->vertices().size(); i++){
      InterpolationElements& elems = _weights[i];
      size_t outOffset = i * dimensions;
      for (query::InterpolationElement& elem : elems) {
        size_t inOffset = (size_t)elem.element->getID() * dimensions;
        for (int dim=0; dim < dimensions; dim++){
          assertion(outOffset + dim < (size_t)outValues.size());
          assertion(inOffset + dim < (size_t)inValues.size());
          outValues(outOffset + dim) += elem.weight * inValues(inOffset + dim);
        }
      }
    }
  }
  else {
    assertion(getConstraint() == CONSERVATIVE, getConstraint());
    DEBUG("Map conservative");
    assertion(_weights.size() == input()->vertices().size(),
               _weights.size(), input()->vertices().size());
    for (size_t i=0; i < input()->vertices().size(); i++){
      size_t inOffset = i * dimensions;
      InterpolationElements& elems = _weights[i];
      for (query::InterpolationElement& elem : elems) {
        size_t outOffset = (size_t)elem.element->getID() * dimensions;
        for ( int dim=0; dim < dimensions; dim++ ){
          assertion(outOffset + dim < (size_t)outValues.size());
          assertion(inOffset + dim < (size_t)inValues.size());
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

  if (getConstraint() == CONSISTENT){
    for(mesh::Vertex& v : input()->vertices()){
      for (size_t i=0; i < output()->vertices().size(); i++) {
        const InterpolationElements& elems = _weights[i];
        for (const query::InterpolationElement& elem : elems) {
          if (elem.element->getID()==v.getID() && elem.weight!=0.0) {
            v.tag();
          }
        }
      }
    }
  }
  else {
    assertion(getConstraint() == CONSERVATIVE, getConstraint());
    for(mesh::Vertex& v : output()->vertices()){
      for (size_t i=0; i < input()->vertices().size(); i++) {
        const InterpolationElements& elems = _weights[i];
        for (const query::InterpolationElement& elem : elems) {
          if (elem.element->getID()==v.getID() && elem.weight!=0.0) {
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

}} // namespace precice, mapping
