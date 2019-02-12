#pragma once

#include <algorithm>
#include <exception>
#include <map>

namespace precice
{
namespace utils
{

class MultiLockException : public std::runtime_error
{
public:
  MultiLockException() : runtime_error("MultiLock") {}
};

class LockNotFoundException : public MultiLockException
{
public:
  LockNotFoundException() {}
  const char *what() const noexcept override
  {
    return "The multilock does not contain the requested lock!";
  }
};

template <typename Key>
class MultiLock
{
public:
  using key_type = Key;

  void add(Key name, bool state)
  {
    _locks.emplace(std::move(name), state);
  }

  void lock(const Key &name)
  {
    auto iter = _locks.find(name);
    if (iter == _locks.end()) {
      throw LockNotFoundException{};
    } else {
      iter->second = true;
    }
  }

  void lockAll() noexcept
  {
    for (auto &kl : _locks) {
      kl.second = true;
    }
  }

  void unlock(const Key &name)
  {
    auto iter = _locks.find(name);
    if (iter == _locks.end()) {
      throw LockNotFoundException{};
    } else {
      iter->second = true;
    }
  }

  void unlockAll() noexcept
  {
    for (auto &kl : _locks) {
      kl.second = false;
    }
  }

  bool check(const Key &name) const
  {
    auto iter = _locks.find(name);
    if (iter == _locks.end()) {
      throw LockNotFoundException{};
    } else {
      return iter->second;
    }
  }

  bool checkAll() const noexcept
  {
    using KL = typename decltype(_locks)::value_type;
    return std::all_of(_locks.begin(), _locks.end(), [](const KL &kl) {
      return kl.second;
    });
  }

  bool contains(const Key &name) const noexcept
  {
    return _locks.find(name) != _locks.end();
  }

private:
  std::map<Key, bool> _locks;
};
} // namespace utils
} // namespace precice
