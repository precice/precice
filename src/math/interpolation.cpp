#include "math/interpolation.hpp"

namespace precice {
namespace math {
namespace interpolation {

std::vector<InterpolationElement> generateInterpolationElements(
    const mesh::Vertex & /*location*/,
    const mesh::Vertex &element)
{
  return {{element, 1.0}};
}

std::vector<InterpolationElement> generateInterpolationElements(
    const mesh::Vertex &location,
    const mesh::Edge &  element)
{
  auto &A = element.vertex(0);
  auto &B = element.vertex(1);

  const auto bcoords = math::barycenter::calcBarycentricCoordsForEdge(
                           A.getCoords(),
                           B.getCoords(),
                           element.getNormal(),
                           location.getCoords())
                           .barycentricCoords;

  std::vector<InterpolationElement> elems;
  elems.emplace_back(A, bcoords(0));
  elems.emplace_back(B, bcoords(1));
  return elems;
}

std::vector<InterpolationElement> generateInterpolationElements(
    const mesh::Vertex &  location,
    const mesh::Triangle &element)
{
  auto &A = element.vertex(0);
  auto &B = element.vertex(1);
  auto &C = element.vertex(2);

  const auto bcoords = math::barycenter::calcBarycentricCoordsForTriangle(
                           A.getCoords(),
                           B.getCoords(),
                           C.getCoords(),
                           element.getNormal(),
                           location.getCoords())
                           .barycentricCoords;

  std::vector<InterpolationElement> elems;
  elems.emplace_back(A, bcoords(0));
  elems.emplace_back(B, bcoords(1));
  elems.emplace_back(C, bcoords(2));
  return elems;
}

} // namespace interpolation
} // namespace math
} // namespace precice