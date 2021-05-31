#pragma once

#include <algorithm>
#include <exception>
#include <map>

namespace precice {
namespace utils {

class MultiLockException : public std::runtime_error {
public:
  MultiLockException()
      : runtime_error("MultiLock") {}
};

class LockNotFoundException : public MultiLockException {
public:
  LockNotFoundException() {}
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
  void lock(const Key &name)
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
  void unlock(const Key &name)
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
  bool check(const Key &name) const
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
  bool contains(const Key &name) const noexcept
  {
    return _locks.find(name) != _locks.end();
  }

private:
  /// The map that keeps track of the locks and their state.
  std::map<Key, bool> _locks;
};
} // namespace utils
} // namespace precice
