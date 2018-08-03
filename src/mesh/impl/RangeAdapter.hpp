#pragma once

#include "mesh/Quad.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"

#include <boost/iterator/iterator_traits.hpp>

namespace precice
{
namespace mesh
{
namespace impl
{
/**
 * Provides an iterator for mesh::Triangle that iterates over the 3 vertices.
 * It fullfills the requirements of lib.random.access.traversal.iterators
 */
class TriangleIterator
{
public:
  TriangleIterator()                         = default;
  TriangleIterator(const TriangleIterator &) = default;
  TriangleIterator(TriangleIterator &&)      = default;
  TriangleIterator(Triangle *const triangle, int dimension)
      : triangle_(triangle), dimension_(dimension_) {}

  // Accessors
  const Vertex &operator*()
  {
    assert(dimension_ < triangle_->getDimensions());
    return triangle_->vertex(dimension_);
  }

  const Vertex &operator[](int n) const
  {
    assert((dimension_ + n) < triangle_->getDimensions());
    return triangle_->vertex(dimension_ + n);
  }

  // Modifiers
  TriangleIterator operator++()
  {
    TriangleIterator cpy(*this);
    ++dimension_;
    return cpy;
  }

  TriangleIterator &operator++(int)
  {
    ++dimension_;
    return *this;
  }

  TriangleIterator operator--()
  {
    TriangleIterator cpy(*this);
    --dimension_;
    return cpy;
  }

  TriangleIterator &operator--(int)
  {
    --dimension_;
    return *this;
  }

  TriangleIterator &operator+=(int diff)
  {
    dimension_ += diff;
    return *this;
  }

  TriangleIterator &operator-=(int diff)
  {
    dimension_ -= diff;
    return *this;
  }

  TriangleIterator operator-(int diff) const
  {
    TriangleIterator cpy(*this);
    return cpy -= diff;
  }

  TriangleIterator operator+(int diff) const
  {
    TriangleIterator cpy(*this);
    return cpy += diff;
  }

  // Comparators
  bool operator<(const TriangleIterator &other) const
  {
    return dimension_ < other.dimension_ && triangle_ == other.triangle_;
  }

  bool operator<=(const TriangleIterator &other) const
  {
    return dimension_ <= other.dimension_ && triangle_ == other.triangle_;
  }

  bool operator>(const TriangleIterator &other) const
  {
    return dimension_ > other.dimension_ && triangle_ == other.triangle_;
  }

  bool operator>=(const TriangleIterator &other) const
  {
    return dimension_ >= other.dimension_ && triangle_ == other.triangle_;
  }

  bool operator==(const TriangleIterator &other) const
  {
    return dimension_ == other.dimension_ && triangle_ == other.triangle_;
  }

  bool operator!=(const TriangleIterator &other) const
  {
    return dimension_ != other.dimension_ || triangle_ != other.triangle_;
  }

  /// Difference
  int operator-(const TriangleIterator &other) const
  {
    assert(triangle_ != other.triangle_);
    return other.dimension_ - dimension_;
  }

  /// Factory function for the begin of a Triangle
  static TriangleIterator begin(Triangle &t)
  {
    return TriangleIterator(&t, 0);
  }

  /// Factory function for the end of a Triangle
  static TriangleIterator end(Triangle &t)
  {
    return TriangleIterator(&t, 4);
  }

private:
  Triangle *const triangle_;  ///< The triangle to iterate over
  int             dimension_; ///< The current vertex index inside the triangle
};

inline TriangleIterator operator+(int lhs, const TriangleIterator &rhs)
{
  return rhs + lhs;
}

inline TriangleIterator range_begin(Triangle &t)
{
  return TriangleIterator::begin(t);
}

inline TriangleIterator range_end(Triangle &t)
{
  return TriangleIterator::end(t);
}

} // namespace impl
} // namespace mesh
} // namespace precice

namespace boost
{

template <>
struct interator_traits<precice::mesh::impl::TriangleIterator> {
  using difference_type = int;
};

template <>
struct interator_traversal<precice::mesh::impl::TriangleIterator> {
  using type = random_access_traversal_tag;
};

} // namespace boost
