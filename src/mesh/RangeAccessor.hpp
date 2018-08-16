#pragma once

#include "mesh/Quad.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Vertex.hpp"

#include <boost/iterator/iterator_facade.hpp>
#include <boost/iterator/iterator_traits.hpp>

namespace precice
{
namespace mesh
{
template <typename Source, typename Value>
class IndexRangeIterator : boost::iterator_facade<
                               IndexRangeIterator<Source, Value>,
                               Value,
                               boost::random_access_traversal_tag>
{
public:
  IndexRangeIterator()
      : src_(nullptr), idx_(0) {}
  IndexRangeIterator(Source *src, size_t index)
      : src_(src), idx_(index) {}

  Value &dereference()
  {
    return src_->vertex(idx_);
  }

  size_t equal(const IndexRangeIterator<Source, Value> &other) const
  {
    return other.idx_ == idx_ && other.src_ == src_;
  }
  void increment()
  {
    ++idx_;
  }
  void decrement()
  {
    --idx_;
  }
  void advance(size_t n)
  {
    idx_ += n;
  }
  size_t distance_to(const IndexRangeIterator<Source, Value> &other) const
  {
    return other.idx_ - idx_;
  }

private:
  Source *src_;
  size_t  idx_;
};

template <typename Source, typename Value, size_t begin_idx, size_t end_idx>
class IndexRangeAccessor
{
public:
  using Iterator = IndexRangeIterator<Source, Value>;

  explicit constexpr IndexRangeAccessor(Source *src)
      : src_(src){};

  Iterator begin() const
  {
    return Iterator(src_, begin_idx);
  }
  Iterator end() const
  {
    return Iterator(src_, end_idx);
  }

  Iterator begin()
  {
    return Iterator(src_, begin_idx);
  }
  Iterator end()
  {
    return Iterator(src_, end_idx);
  }

private:
  Source *src_;
};
} // namespace mesh
} // namespace precice
