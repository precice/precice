#pragma once

#include <boost/test/unit_test.hpp>
#include <numeric>
#include <stdexcept>
#include "math/math.hpp"
#include "precice/config/Configuration.hpp"
#include "utils/ManageUniqueIDs.hpp"
#include "utils/Parallel.hpp"
// #include "utils/traits.hpp"
#include "xml/XMLTag.hpp"

namespace precice {
namespace testing {

namespace bt = boost::unit_test;
using Par    = precice::utils::Parallel;

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
    Participants participants{
        "Serial"_on(ranks)};
    initialize(participants);
  }

  template <class... T>
  TestContext(Ranks ranks, T... args)
      : _simple(true)
  {
    Participants participants{
        "Serial"_on(ranks)};
    //static_assert(
    //    precice::utils::conjunction< std::is_base_of<OptionTag, T>::value...>::value,
    //    "Remaining arguments to a serial test case have to be options");
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

  m2n::PtrM2N connect(const std::string &acceptor, const std::string &requestor) const;

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

namespace inject {
using precice::testing::Require;
using precice::testing::operator""_rank;
using precice::testing::operator""_ranks;
using precice::testing::operator""_on;
} // namespace inject

#define PRECICE_TEST(...)                             \
  using namespace precice::testing::inject;           \
  precice::testing::TestContext context{__VA_ARGS__}; \
  if (context.invalid) {                              \
    return;                                           \
  }

/// Boost.Test decorator that unconditionally deletes the test.
class Deleted : public bt::decorator::base {
private:
  virtual void apply(bt::test_unit &tu)
  {
    bt::framework::get<bt::test_suite>(tu.p_parent_id).remove(tu.p_id);
  }

  virtual bt::decorator::base_ptr clone() const
  {
    return bt::decorator::base_ptr(new Deleted());
  }
};

/// struct giving access to the impl of a befriended class or struct
struct WhiteboxAccessor {
  /** Returns the impl of the obj by reference.
     *
     * Returns a reference to the object pointed to by the _impl of a class.
     * This class needs to be friend of T.
     *
     * @param[in] obj The object to fetch the impl from.
     * @returns a lvalue reference to the impl object.
     */
  template <typename T>
  static auto impl(T &obj) -> typename std::add_lvalue_reference<decltype(*(obj._impl))>::type
  {
    return *obj._impl;
  }
};

/// equals to be used in tests. Prints both operatorans on failure
template <class DerivedA, class DerivedB>
boost::test_tools::predicate_result equals(const Eigen::MatrixBase<DerivedA> &A,
                                           const Eigen::MatrixBase<DerivedB> &B,
                                           double                             tolerance = math::NUMERICAL_ZERO_DIFFERENCE)
{
  if (not math::equals(A, B, tolerance)) {
    boost::test_tools::predicate_result res(false);
    Eigen::IOFormat                     format;
    if (A.cols() == 1) {
      format.rowSeparator = ", ";
      format.rowPrefix    = "  ";
    }
    res.message() << "\n"
                  << A.format(format) << " != \n"
                  << B.format(format);
    return res;
  }
  return true;
}

/// equals to be used in tests. Prints both operatorans on failure
template <class Scalar>
typename std::enable_if<std::is_arithmetic<Scalar>::value, boost::test_tools::predicate_result>::type equals(const Scalar a, const Scalar b, const Scalar tolerance = math::NUMERICAL_ZERO_DIFFERENCE)
{
  if (not math::equals(a, b, tolerance)) {
    boost::test_tools::predicate_result res(false);
    res.message() << "Not equal: " << a << "!=" << b;
    return res;
  }
  return true;
}

/// Returns $PRECICE_ROOT/src, the base path to the sources.
std::string getPathToSources();

/** Generates a new mesh id for use in tests.
 *
 * @returns a new unique mesh ID
 */
inline int nextMeshID()
{
  static utils::ManageUniqueIDs manager;
  return manager.getFreeID();
}

} // namespace testing
} // namespace precice
