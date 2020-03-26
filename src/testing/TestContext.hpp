#pragma once

#include <string>
#include <vector>

#include "m2n/M2N.hpp"
#include "utils/Parallel.hpp"

namespace precice {
namespace testing {

struct Ranks {
  int value;
};

inline constexpr Ranks operator""_ranks(unsigned long long value)
{
  return Ranks{static_cast<int>(value)};
}

inline constexpr Ranks operator""_rank(unsigned long long value)
{
  return (value == 1) ? Ranks{1} : throw;
}

struct Participant {
  std::string name;
  int         size   = 1;
  bool        initMS = false;

  explicit Participant(std::string n)
      : name(std::move(n)){};

  Participant &operator()(Ranks rsize)
  {
    size = rsize.value;
    return *this;
  }

  Participant &setupMasterSlaves()
  {
    initMS = true;
    return *this;
  }
};

inline Participant operator""_on(const char *name, long unsigned int)
{
  return Participant{name};
}

static_assert(std::is_same<Participant &, decltype(""_on(2_ranks))>::value, "");

enum class Require {
  PETSc,
  Events,
};

enum struct ConnectionType {
  GatherScatter,
  PointToPoint
};

struct ConnectionOptions {
  ConnectionOptions()             = default;
  bool           useOnlyMasterCom = false;
  bool           useTwoLevelInit  = false;
  ConnectionType type             = ConnectionType::GatherScatter;
};

class TestContext {
public:
  using Participants = std::vector<Participant>;

  std::string name;
  int         rank    = 0;
  int         size    = 1;
  bool        invalid = false;

  TestContext() = default;

  template <class... T>
  TestContext(Ranks ranks)
      : _simple(true)
  {
    Participants participants{"Serial"_on(ranks)};
    initialize(participants);
  }

  template <class... T>
  TestContext(Ranks ranks, T... args)
      : _simple(true)
  {
    Participants participants{"Serial"_on(ranks)};
    handleOptions(participants, args...);
    initialize(participants);
  }

  template <class... T>
  TestContext(T... args)
  {
    Participants participants;
    handleOptions(participants, args...);
    initialize(participants);
  }

  ~TestContext() noexcept;

  bool hasSize(int size) const;

  bool isNamed(const std::string &name) const;

  bool isRank(int rank) const;

  bool isMaster() const;

  m2n::PtrM2N connect(const std::string &acceptor, const std::string &requestor, const ConnectionOptions &options = ConnectionOptions{}) const;

  std::string describe() const;

private:
  bool _petsc  = false;
  bool _events = false;
  bool _simple = false;
  bool _initMS = false;

  utils::Parallel::CommStatePtr _contextComm;

  std::vector<std::string> _names;

  void handleOption(Participants &participants, Participant participant);
  void handleOption(Participants &participants, testing::Require requirement);

  template <class LastOption>
  void handleOptions(Participants &participants, LastOption &last)
  {
    handleOption(participants, last);
  }

  template <class NextOption, class... Rest>
  void handleOptions(Participants &participants, NextOption &next, Rest &... rest)
  {
    handleOption(participants, next);
    handleOptions(participants, rest...);
  }

  void setContextFrom(const Participant &p, int rank);

  void initialize(const Participants &participants);
  void initializeMPI(const Participants &participants);
  void initializeMasterSlave();
  void initializePetsc();
  void initializeEvents();
};

} // namespace testing
} // namespace precice
