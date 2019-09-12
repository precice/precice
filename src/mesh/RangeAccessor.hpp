#pragma once

#include <boost/iterator/iterator_facade.hpp>

namespace precice {
namespace mesh {
/** random-access iterator over an indexable Source.
 * 
 * @tparam Source the underlying container to index into
 * @tparam Value the resulting value
 *
 * @note This version currently only supports Sources with a const `src.vertex(index).getCoords()` access.
 */
template <typename Source, typename Value>
class IndexRangeIterator : public boost::iterator_facade<
                               IndexRangeIterator<Source, const Value>,
                               const Value,
                               boost::random_access_traversal_tag> {
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
  Source *src_{nullptr}; ///< the source to access
  size_t  idx_{0};       ///< the current index to access
};

} // namespace mesh
} // namespace precice
