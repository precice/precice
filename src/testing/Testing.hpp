#pragma once

#include <Eigen/Core>
#include <boost/test/unit_test.hpp>
#include <string>
#include <type_traits>

#include "math/differences.hpp"
#include "math/math.hpp"
#include "testing/TestContext.hpp"
#include "utils/IntraComm.hpp"

namespace precice::testing {

namespace bt = boost::unit_test;

DataID operator"" _dataID(unsigned long long n);

namespace inject {
using precice::testing::Require;
using precice::testing::operator""_rank;
using precice::testing::operator""_ranks;
using precice::testing::operator""_on;
using precice::testing::operator""_dataID;
} // namespace inject

#define PRECICE_TEST(...)                             \
  using namespace precice::testing::inject;           \
  precice::testing::TestContext context{__VA_ARGS__}; \
  if (context.invalid) {                              \
    return;                                           \
  }                                                   \
  BOOST_TEST_MESSAGE(context.describe());             \
  boost::unit_test::framework::add_context(BOOST_TEST_LAZY_MSG(context.describe()), true);

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

/// equals to be used in tests. Compares two matrices using a given tolerance. Prints both operands on failure.
template <class DerivedA, class DerivedB>
boost::test_tools::predicate_result equals(const Eigen::MatrixBase<DerivedA> &A,
                                           const Eigen::MatrixBase<DerivedB> &B,
                                           double                             tolerance = math::NUMERICAL_ZERO_DIFFERENCE)
{
  auto approx = [tolerance](double a, double b) -> bool {
    if (std::max(std::abs(a), std::abs(b)) < tolerance) {
      return true;
    } else {
      return std::abs(a - b) <= tolerance * std::max(std::abs(a), std::abs(b));
    }
  };

  if (!A.binaryExpr(B, approx).all()) {
    boost::test_tools::predicate_result res(false);
    Eigen::IOFormat                     format;
    if (A.cols() == 1) {
      format.rowSeparator = ", ";
      format.rowPrefix    = "  ";
      format.precision    = Eigen::FullPrecision;
    }
    res.message() << "\n"
                  << A.format(format) << " != \n"
                  << B.format(format) << '\n'
                  << (A - B).cwiseAbs().format(format) << " difference\n"
                  << A.binaryExpr(B, approx).format(format) << "in tolerance (" << tolerance << ')';
    return res;
  }
  return true;
}

/// equals to be used in tests. Compares two std::vectors using a given tolerance. Prints both operands on failure
boost::test_tools::predicate_result equals(const std::vector<float> &VectorA,
                                           const std::vector<float> &VectorB,
                                           float                     tolerance = math::NUMERICAL_ZERO_DIFFERENCE);

boost::test_tools::predicate_result equals(const std::vector<double> &VectorA,
                                           const std::vector<double> &VectorB,
                                           double                     tolerance = math::NUMERICAL_ZERO_DIFFERENCE);

/// equals to be used in tests. Compares two scalar numbers using a given tolerance. Prints both operands on failure
boost::test_tools::predicate_result equals(float a, float b, float tolerance = math::NUMERICAL_ZERO_DIFFERENCE);

boost::test_tools::predicate_result equals(double a, double b, double tolerance = math::NUMERICAL_ZERO_DIFFERENCE);

/// Returns the base path of the repo.
std::string getPathToRepository();

/// Returns the base path to the sources.
std::string getPathToSources();

/// Returns the base path to the integration tests.
std::string getPathToTests();

/// Returns the name of the current test.
std::string getTestName();

/// Returns the full path to the file containing the current test.
std::string getTestPath();

/** Generates a new mesh id for use in tests.
 *
 * @returns a new unique mesh ID
 */
int nextMeshID();

} // namespace precice::testing
