#pragma once

#include <boost/iterator/indirect_iterator.hpp>
#include <vector>

#include "utils/assertion.hpp"

namespace precice {
namespace utils {

/// Wrapper around std::vector for transparent handling of pointer objects.
template <typename CONTENT_T>
class ptr_vector {
public:
  /// Type of wrapped vector container.
  typedef std::vector<CONTENT_T *> container;
  typedef CONTENT_T                value_type; // necessary to be standard conform

  // The size_type of the wrapped vector
  using size_type = typename container::size_type;

  /// Type of iterator hiding pointers.
  typedef boost::indirect_iterator<
      typename container::iterator,
      CONTENT_T *,
      boost::use_default,
      CONTENT_T &>
      iterator;

  /// Type of const_iterator hiding pointers.
  typedef boost::indirect_iterator<
      typename container::const_iterator,
      CONTENT_T *,
      boost::use_default,
      CONTENT_T &>
      const_iterator;

  /// Returns iterator to first element in vector.
  iterator begin()
  {
    return iterator(_content.begin());
  }

  /**
    * @brief Returns const interator to first element in vector.
    */
  const_iterator begin() const
  {
    return const_iterator(_content.begin());
  }

  /**
    * @brief Returns iterator behind last element in vector.
    */
  iterator end()
  {
    return iterator(_content.end());
  }

  /**
    * @brief Returns const_iterator behind last element in vector.
    */
  const_iterator end() const
  {
    return const_iterator(_content.end());
  }

  /**
    * @brief Returns count of elements in vector.
    */
  size_t size() const
  {
    return _content.size();
  }

  /**
    * @brief Returns reference to element with given index [0, count[.
    */
  CONTENT_T &operator[](size_t index)
  {
    PRECICE_ASSERT(index < _content.size());
    return *_content[index];
  }

  /**
    * @brief Returns const reference to element with given index [0, count[.
    */
  const CONTENT_T &operator[](size_t index) const
  {
    PRECICE_ASSERT(index < _content.size());
    return *_content[index];
  }

  CONTENT_T &back()
  {
    return *_content.back();
  }

  const CONTENT_T &back() const
  {
    return *_content.back();
  }

  /**
    * @brief Adds element to the end of the vector.
    */
  void push_back(CONTENT_T *content)
  {
    _content.push_back(content);
  }

  /**
    * @brief Inserts elements into vector.
    *
    * @param position [IN] Iterator at position before which the elements are
    *        inserted.
    * @param from [IN] Iterator at position of first element to be inserted.
    * @param to [IN] Iterator at position behind last element to be inserted.
    */
  void insert(iterator position, iterator from, iterator to)
  {
    _content.insert(position.base(), from.base(), to.base());
  }

  /**
    * @brief Removes all pointers to elements from the vector (no deletetion).
    */
  void clear()
  {
    _content.clear();
  }

  /**
    * @brief Returns true, if no pointers are contained in the vector.
    */
  bool empty() const
  {
    return _content.empty();
  }

  void deleteElements()
  {
    for (CONTENT_T *elem : _content) {
      PRECICE_ASSERT(elem != NULL);
      delete (elem);
    }
  }

private:
  container _content;
};

} // namespace utils
} // namespace precice
