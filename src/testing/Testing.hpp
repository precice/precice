#pragma once

#include <Eigen/Core>
#include <boost/test/unit_test.hpp>
#include <string>
#include <type_traits>
#include "math/differences.hpp"
#include "math/math.hpp"
#include "testing/TestContext.hpp"
#include "utils/ManageUniqueIDs.hpp"
#include "utils/MasterSlave.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {
namespace testing {

namespace bt = boost::unit_test;

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
  }                                                   \
  BOOST_TEST_MESSAGE(context.describe());

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

/// Holds rank, owner, position, value of a single vertex
struct VertexSpecification {
  int                 rank;
  int                 owner;
  std::vector<double> position;
  std::vector<double> value;

  const Eigen::Map<const Eigen::VectorXd> asEigen() const
  {
    return Eigen::Map<const Eigen::VectorXd>{position.data(), static_cast<Eigen::Index>(position.size())};
  }
};

/// Holds the indices of the vertices of the edge and the rank in which the edge is going to be created
struct EdgeSpecification {
  std::vector<int> vertices;
  int              rank;
};

/// Holds the indices of the edges of the face and the rank in which the face is going to be created
struct FaceSpecification {
  std::vector<int>              edges;
  int                           rank;
};

/*
- -1 on rank means all ranks
- -1 on owner rank means no rank
- x, y, z is position of vertex, z is optional, 2D mesh will be created then
*/
struct MeshSpecification {
  std::vector<VertexSpecification> vertices;
  std::vector<EdgeSpecification>   edges;
  std::vector<FaceSpecification>   faces;

  MeshSpecification(std::vector<VertexSpecification> vertices_)
      : vertices(vertices_) {}

  MeshSpecification(std::vector<VertexSpecification> vertices_, std::vector<EdgeSpecification> edges_)
      : vertices(vertices_), edges(edges_) {}

  MeshSpecification(std::vector<VertexSpecification> vertices_, std::vector<EdgeSpecification> edges_, std::vector<FaceSpecification> faces_)
      : vertices(vertices_), edges(edges_), faces(faces_) {}
};

/*
ReferenceSpecification format:
{ {rank, {v}, ... }
- -1 on rank means all ranks
- v is the expected value of n-th vertex on that particular rank
*/
using ReferenceSpecification = std::vector<std::pair<int, std::vector<double>>>;

} // namespace testing
} // namespace precice
