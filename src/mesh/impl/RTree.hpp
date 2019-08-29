#pragma once

#include "mesh/RTree.hpp"

namespace precice {
namespace mesh {
namespace impl {


/** Holding a reference to a Mesh it is a functor for packing
 *  PrimitiveIndex into AABBs.
 */
class AABBGenerator
{
public:
  /// Constructs a generator for a given mesh
  explicit AABBGenerator(const mesh::Mesh &mesh)
      : mesh_(mesh){};

  /** constructs a AABB of a PrimitiveIndex.
   *
   * This operator takes a PrimitiveIndex and fetches the actual primitive from the mesh.
   * It then constructs an AABB using boost::geometry::return_envelope() and returns it.
   *
   * @param pi the index to build the AABB for
   * @returns a AABB for the primitive
   */
  AABB operator()(const PrimitiveIndex &pi) const
  {
    switch (pi.type) {
    case (Primitive::Vertex):
      return boost::geometry::return_envelope<AABB>(mesh_.vertices()[pi.index]);
    case (Primitive::Edge):
      return boost::geometry::return_envelope<AABB>(mesh_.edges()[pi.index]);
    case (Primitive::Triangle):
      return boost::geometry::return_envelope<AABB>(mesh_.triangles()[pi.index]);
    case (Primitive::Quad):
      return boost::geometry::return_envelope<AABB>(mesh_.quads()[pi.index]);
    default:
      throw std::invalid_argument{"AABB generation for this Primitive is not implemented!"};
    }
  }

private:
  /// The mesh used look the primitives up
  const mesh::Mesh &mesh_;
};

/** indexes a container of primitives into an rtree.
 *
 * The algorithm indexes every primitive of the given container using the passed generator.
 * It inserts its result and the PrimitiveIndex into the rtree as a std::pair.
 *
 * @param[IN, OUT] rtree the rtree to index into
 * @param gen the Generator generating something to index rtree::value_type.
 * @param conti the Container to index
 */
template <typename Container, typename Generator = AABBGenerator>
void indexPrimitive(PrimitiveRTree &rtree, const Generator& gen, const Container &conti)
{
  using ValueType = typename std::remove_reference<typename std::remove_cv<typename Container::value_type>::type>::type;
  for (size_t i = 0; i < conti.size(); ++i) {
    PrimitiveIndex index{as_primitive_enum<ValueType>::value, i};
    rtree.insert(std::make_pair(gen(index), index));
  }
}

}}}
