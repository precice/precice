#pragma once

#include <boost/iterator/iterator_facade.hpp>

namespace precice
{
namespace mesh
{
template <typename Source, typename Value>
class IndexRangeIterator : public boost::iterator_facade<
                               IndexRangeIterator<Source, const Value>,
                               const Value,
                               boost::random_access_traversal_tag>
{
public:
  IndexRangeIterator() = default;
  IndexRangeIterator(Source *src, size_t index)
      : src_(src), idx_(index) {}

  const Value &dereference() const
  {
    using Coord = decltype(src_->vertex(idx_).getCoords());
    static_assert(
            std::is_reference<Coord>::value,
            "Coordinate type must be a reference!");
    static_assert(
            std::is_convertible<Coord, Value>::value,
            "Exposed and accessed types must match!");
    return static_cast<const Value &>(src_->vertex(idx_).getCoords());
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
  Source *src_{nullptr};
  size_t  idx_{0};
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
