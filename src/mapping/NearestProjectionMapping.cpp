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

void NearestProjectionMapping:: computeMapping()
{
  TRACE()
  precice::utils::Event e("map.np.computeMapping.From" + input()->getName() + "To" + output()->getName(), precice::syncMode);
  if (getConstraint() == CONSISTENT){
      computeMappingConsistent();
  } else {
      computeMappingConservative();
  }
  auto negativeWeights = std::count_if(_weights.begin(), _weights.end(), [](const InterpolationElements& elems) -> bool {
          return std::any_of(elems.begin(), elems.end(), [](const query::InterpolationElement& elem) -> bool {
                  return elem.weight < 0.0;
                  });
          });
  INFO("Extrapolating on " << negativeWeights << " / " << _weights.size() << " associaltions");
}

void NearestProjectionMapping:: computeMappingConsistent()
{
    TRACE(input()->vertices().size(), output()->vertices().size());
    DEBUG("Compute consistent mapping");
    assertion(getConstraint() == CONSISTENT, getConstraint());
    const auto &oVertices = output()->vertices();
    _weights.resize(oVertices.size());

    if (getDimensions() == 2) {
        const auto &iEdges    = input()->edges();
        CHECK(oVertices.empty() || !iEdges.empty(), "Mesh \"" << input()->getName() << "\" does not contain Edges to project onto.");

        auto rtree = mesh::rtree::getEdgeRTree(input());
        for (size_t i = 0; i < oVertices.size(); i++) {
            const Eigen::VectorXd &coords = oVertices[i].getCoords();
            // Search for the output vertex inside the input mesh
            rtree->query(boost::geometry::index::nearest(coords, 1),
                     boost::make_function_output_iterator([&](size_t eid) {
                         auto& weights = _weights[i];
                         weights = query::generateInterpolationElements(oVertices[i], iEdges[eid]);
                         CHECK(!weights.empty(), "No interpolation elements for current vertex!");
                         if (std::any_of(weights.begin(), weights.end(), [](const query::InterpolationElement& elem){return elem.weight < 0;})) {
                            WARN("Mapping \"" << output()->getName() << "\" contains vertex (" << oVertices[i] << ") which has negative weights indicating non-matching meshes!");
                         }
                    }));
        }
    } else {
        const auto &iTriangles = input()->triangles();
        CHECK(oVertices.empty() || !iTriangles.empty(), "Mesh \"" << input()->getName() << "\" does not contain Triangles to project onto.");

        mesh::rtree::clear(*input());
        auto rtree = mesh::rtree::getTriangleRTree(input());
        DEBUG("RTree contains " << rtree->size() << " elements");
        CHECK(rtree->size() == iTriangles.size(), "Caching is wrong");
        using IndexType = typename mesh::rtree::triangle_traits::IndexType;
        using MatchType = std::pair<double, int>;
        std::vector<MatchType> matches;
        matches.reserve(10);
        for (size_t i = 0; i < oVertices.size(); i++) {
            const Eigen::VectorXd &coords = oVertices[i].getCoords();
            // Search for the output vertex inside the input mesh
            matches.clear();
            rtree->query(boost::geometry::index::nearest(coords, 10),
                    boost::make_function_output_iterator([&](const IndexType& match) {
                        using boost::geometry::distance;
                        auto dist = distance(coords, iTriangles[match.second]);
                        matches.emplace_back(dist, match.second);
                    }));
            DEBUG("Matches for " << i << " " << oVertices[i] << " :");
            std::sort(matches.begin(), matches.end(), [&](const MatchType& lhs, const MatchType& rhs) { return lhs.first < rhs.first; });
            for(auto& match: matches) {
                DEBUG("match with dist " << match.first << " is "  << match.second << ":" << iTriangles[match.second]);
            }
            auto closest = *std::min_element(matches.begin(), matches.end(), 
                    [&](const MatchType& lhs, const MatchType& rhs) {
                        return lhs.first < rhs.first;
                    });
             DEBUG("closest match with dist " << closest.first << " : " << closest.second << ":" << iTriangles[closest.second]);
             auto& weights = _weights[i];
             weights = query::generateInterpolationElements(oVertices[i], iTriangles[closest.second]);
             DEBUG("weights: " << weights);
             CHECK(!weights.empty(), "No interpolation elements for current vertex!");
             if (std::any_of(weights.begin(), weights.end(), [](const query::InterpolationElement& elem){return elem.weight < 0;})) {
                WARN("Mapping \"" << output()->getName() << "\" contains vertex (" << oVertices[i] << ") which has negative weights indicating non-matching meshes!");
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
    const auto &iVertices = input()->vertices();
    _weights.resize(iVertices.size());

    if (getDimensions() == 2) {
        const auto &oEdges    = output()->edges();
        CHECK(iVertices.empty() || !oEdges.empty(), "Mesh \"" << output()->getName() << "\" does not contain Edges to project onto.");

        auto rtree = mesh::rtree::getEdgeRTree(output());
        for (size_t i = 0; i < iVertices.size(); i++) {
            const Eigen::VectorXd &coords = iVertices[i].getCoords();
            // Search for the output vertex inside the input mesh
            rtree->query(boost::geometry::index::nearest(coords, 1),
                     boost::make_function_output_iterator([&](size_t eid) {
                         auto& weights = _weights[i];
                         weights = query::generateInterpolationElements(iVertices[i], oEdges[eid]);
                         DEBUG("match for " << iVertices[i] << " is " << oEdges[eid] << " with weight " << weights);
                         CHECK(!weights.empty(), "No interpolation elements for current vertex!");
                         if (std::any_of(weights.begin(), weights.end(), [](const query::InterpolationElement& elem){return elem.weight < 0;})) {
                            WARN("Mesh \"" << input()->getName() << "\" contains vertex (" << iVertices[i] << ") which has negative weights indicating non-matching meshes!");
                         }
                    }));
        }
    } else {
        const auto &oTriangles = output()->triangles();
        CHECK(iVertices.empty() || !oTriangles.empty(), "Mesh \"" << output()->getName() << "\" does not contain Triangles to project onto.");

        auto rtree = mesh::rtree::getTriangleRTree(output());
        using IndexType = typename mesh::rtree::triangle_traits::IndexType;
        using MatchType = std::pair<double, int>;
        std::vector<MatchType> matches;
        for (size_t i = 0; i < iVertices.size(); i++) {
            const Eigen::VectorXd &coords = iVertices[i].getCoords();
            // Search for the output vertex inside the input mesh
            rtree->query(boost::geometry::index::nearest(coords, 10),
                     boost::make_function_output_iterator([&](const IndexType& match) {
                    using boost::geometry::distance;
                    auto dist = distance(coords, oTriangles[match.second]);
                    matches.emplace_back(dist, match.second);
                    }));
            DEBUG("Matches for " << iVertices[i] << " :");
            for(auto& match: matches) {
                DEBUG("match with dist " << match.first << " : " << oTriangles[match.second]);
            }
            auto closest = *std::min_element(matches.begin(), matches.end(), 
                    [&](const MatchType& lhs, const MatchType& rhs) {
                        return lhs.first < rhs.first;
                    });
             DEBUG("closest match with dist " << closest.first << " : " << oTriangles[closest.second]);
             auto& weights = _weights[i];
             weights = query::generateInterpolationElements(iVertices[i], oTriangles[closest.second]);
             DEBUG("weights: " << weights);
             CHECK(!weights.empty(), "No interpolation elements for current vertex!");
             if (std::any_of(weights.begin(), weights.end(), [](const query::InterpolationElement& elem){return elem.weight < 0;})) {
                 WARN("Mesh \"" << input()->getName() << "\" contains vertex (" << iVertices[i] << ") which has negative weights indicating non-matching meshes!");
             }
             matches.clear();
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
  DEBUG("Compute Mapping for Tagging");

  computeMapping();
  DEBUG("Tagging First Round");

  std::set<int> tagged;

  if (getConstraint() == CONSISTENT){
    for(mesh::Vertex& v : input()->vertices()){
      for (size_t i=0; i < output()->vertices().size(); i++) {
        const InterpolationElements& elems = _weights[i];
        for (const query::InterpolationElement& elem : elems) {
          if (elem.element->getID()==v.getID() && !math::equals(elem.weight,0.0)) {
            v.tag();
            tagged.insert(v.getID());
          }
        }
      }
    }
    DEBUG("First Round Tagged " << tagged.size() << "/" << input()->vertices().size() << " Vertices");
  }
  else {
    assertion(getConstraint() == CONSERVATIVE, getConstraint());
    for(mesh::Vertex& v : output()->vertices()){
      for (size_t i=0; i < input()->vertices().size(); i++) {
        const InterpolationElements& elems = _weights[i];
        for (const query::InterpolationElement& elem : elems) {
          if (elem.element->getID()==v.getID() && !math::equals(elem.weight,0.0)) {
            v.tag();
            tagged.insert(v.getID());
          }
        }
      }
    }
    DEBUG("First Round Tagged " << tagged.size() << "/" << output()->vertices().size() << " Vertices");
  }

  clear();
}

void NearestProjectionMapping::tagMeshSecondRound()
{
  TRACE();
  // for NP mapping no operation needed here
}

}} // namespace precice, mapping
