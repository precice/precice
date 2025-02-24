#pragma once

#include <algorithm>
#include <exception>
#include <map>

namespace precice::utils {

class MultiLockException : public std::runtime_error {
public:
  MultiLockException()
      : runtime_error("MultiLock") {}
};

class LockNotFoundException : public MultiLockException {
public:
  LockNotFoundException() = default;
  const char *what() const noexcept override
  {
    return "The multilock does not contain the requested lock!";
  }
};

/// Class handling multiple locks allowing global lock and unlock operations.
template <typename Key>
class MultiLock {
public:
  /// The type of the key
  using key_type = Key;

  /// The size type
  using size_type = std::size_t;

  /** @brief Adds a lock with a given state
   *
   * Adding an already existent lock does nothing.
   *
   * @param[in] name the name of the lock
   * @param[in] state the initial state of the lock
   */
  void add(Key name, bool state)
  {
    _locks.emplace(std::move(name), state);
  }

  /** @brief Locks a given lock.
   *
   * @param[in] name the name of to lock
   */
  template <typename K>
  void lock(const K &name)
  {
    auto iter = _locks.find(name);
    if (iter == _locks.end()) {
      throw LockNotFoundException{};
    } else {
      iter->second = true;
    }
  }

  /// Locks all known locks.
  void lockAll() noexcept
  {
    for (auto &kl : _locks) {
      kl.second = true;
    }
  }

  /** @brief Unlocks a given lock.
   *
   * @param[in] name the name of to unlock
   */
  template <typename K>
  void unlock(const K &name)
  {
    auto iter = _locks.find(name);
    if (iter == _locks.end()) {
      throw LockNotFoundException{};
    } else {
      iter->second = false;
    }
  }

  /// Unlocks all known locks.
  void unlockAll() noexcept
  {
    for (auto &kl : _locks) {
      kl.second = false;
    }
  }

  /// Removes all known locks
  void clear() noexcept
  {
    _locks.clear();
  }

  /** @brief Checks the status of a lock
   *
   * @param[in] name the name of the lock to check
   *
   * @returns whether the lock is locked
   */
  template <typename K>
  bool check(const K &name) const
  {
    auto iter = _locks.find(name);
    if (iter == _locks.end()) {
      throw LockNotFoundException{};
    } else {
      return iter->second;
    }
  }

  /// Checks whether all locks are locked.
  bool checkAll() const noexcept
  {
    using KL = typename decltype(_locks)::value_type;
    return std::all_of(_locks.begin(), _locks.end(), [](const KL &kl) {
      return kl.second;
    });
  }

  /** @brief Checks whether a lock is known.
   *
   * @param[in] name the name to check
   *
   * @returns whether the name is a known lock
   */
  template <typename K>
  bool contains(const K &name) const noexcept
  {
    return _locks.find(name) != _locks.end();
  }

  /// Returns the total count of locks
  size_type size() const noexcept
  {
    return _locks.size();
  }

  /// Returns the count of locked locks
  size_t countLocked() const
  {
    return std::count_if(_locks.begin(), _locks.end(), [](typename map_type::value_type const &kv) { return kv.second; });
  }

  /// Returns the count of unlocked locks
  size_t countUnlocked() const
  {
    return std::count_if(_locks.begin(), _locks.end(), [](typename map_type::value_type const &kv) { return not kv.second; });
  }

private:
  using map_type = typename std::map<Key, bool, std::less<>>;

  /// The map that keeps track of the locks and their state.
  map_type _locks;
};
} // namespace precice::utils
