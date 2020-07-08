#pragma once

#include <boost/container/flat_set.hpp>
#include <functional>

namespace precice {
namespace utils {

/// Manages a set of unique IDs.
class ManageUniqueIDs {
public:
  // Returns the next free, i.e. unique, ID.
  int getFreeID();

  /**
    * @brief Inserts an ID which has to be unique.
    *
    * The inserted ID has to be different to all other IDs inserted and obtained
    * from getFreeID().
    */
  bool insertID(int id);

  /// Resets all retrieved and inserted IDs.
  void resetIDs();

private:
  /// Stores all used IDs.
  boost::container::flat_set<int> _ids;

  /// Marks next ID to be given, from lower to higher values.
  int _lowerLimit = 0;
};

} // namespace utils
} // namespace precice
