#include "testing/QuickTest.hpp"
#include "testing/Testing.hpp"

#include <precice/SolverInterface.hpp>

BOOST_AUTO_TEST_SUITE(TestingTests) // Use name of the module, e.g. subdirectory below src/, suffixed with Tests

BOOST_AUTO_TEST_SUITE(Examples) // If your file contains multiple tests, put them in a test suite

BOOST_AUTO_TEST_SUITE(QuickTests)

BOOST_AUTO_TEST_CASE(CheckAtTheEnd)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  using namespace precice::testing;
  std::string config = getPathToRepository() + "/src/testing/tests/quicktest.xml";

  // First initialize the SolverInterface as usual
  precice::SolverInterface interface(context.name, config, context.rank, context.size);

  if (context.isNamed("SolverOne")) {
    // Create a quicktest instance for MeshOne
    QuickTest qt{interface, "MeshOne"};
    qt.setVertices({0.0, 0.0, 0.0,
                    1.0, 0.0, 0.0})
        .initialize()
        .initializeData()
        .read("DataTwo"_vector)
        .write("DataOne"_scalar, {1, 2})
        .advance()
        .read("DataTwo"_vector)
        .finalize();

    auto expected0 = {1, 1, 1, 4, 4, 4};
    auto expected1 = {0, 0, 0, 3, 3, 3};
    BOOST_TEST(qt.reads.at(0) == expected0, boost::test_tools::per_element());
    BOOST_TEST(qt.reads.at(1) == expected1, boost::test_tools::per_element());
  } else {
    //  Create a quicktest instance for MeshTwo
    QuickTest qt{interface, "MeshTwo"};
    qt.setVertices({0.0, 0.0, 0.0,
                    0.33, 0.0, 0.0,
                    0.66, 0.0, 0.0,
                    1.0, 0.0, 0.0})
        .initialize()
        .write("DataTwo"_vector, {1, 1, 1,
                                  2, 2, 2,
                                  3, 3, 3,
                                  4, 4, 4})
        .initializeData()
        .read("DataOne"_scalar)
        .write("DataTwo"_vector, {0, 0, 0,
                                  1, 1, 1,
                                  2, 2, 2,
                                  3, 3, 3})
        .advance()
        .finalize();

    auto expected0 = {1, 1, 2, 2};
    BOOST_TEST(qt.reads.at(0) == expected0, boost::test_tools::per_element());
  }
}

BOOST_AUTO_TEST_CASE(CheckAfterRead)
{
  PRECICE_TEST("SolverOne"_on(1_rank), "SolverTwo"_on(1_rank));
  using namespace precice::testing;
  std::string config = getPathToRepository() + "/src/testing/tests/quicktest.xml";

  // First initialize the SolverInterface as usual
  precice::SolverInterface interface(context.name, config, context.rank, context.size);

  if (context.isNamed("SolverOne")) {
    // Create a quicktest instance for MeshOne
    QuickTest qt{interface, "MeshOne"};
    qt.setVertices({0.0, 0.0, 0.0,
                    1.0, 0.0, 0.0})
        .initialize()
        .initializeData()
        .read("DataTwo"_vector);

    auto expected0 = {1, 1, 1, 4, 4, 4};
    BOOST_TEST(qt.last() == expected0, boost::test_tools::per_element());

    qt.write("DataOne"_scalar, {1, 2})
        .advance()
        .read("DataTwo"_vector);

    auto expected1 = {0, 0, 0, 3, 3, 3};
    BOOST_TEST(qt.last() == expected1, boost::test_tools::per_element());
    qt.finalize();
  } else {
    //  Create a quicktest instance for MeshTwo
    QuickTest qt{interface, "MeshTwo"};
    qt.setVertices({0.0, 0.0, 0.0,
                    0.33, 0.0, 0.0,
                    0.66, 0.0, 0.0,
                    1.0, 0.0, 0.0})
        .initialize()
        .write("DataTwo"_vector, {1, 1, 1,
                                  2, 2, 2,
                                  3, 3, 3,
                                  4, 4, 4})
        .initializeData()
        .read("DataOne"_scalar);
    auto expected0 = {1, 1, 2, 2};
    BOOST_TEST(qt.last() == expected0, boost::test_tools::per_element());
    qt.write("DataTwo"_vector, {0, 0, 0,
                                1, 1, 1,
                                2, 2, 2,
                                3, 3, 3})
        .advance()
        .finalize();
  }
}

BOOST_AUTO_TEST_SUITE_END() // QuickTests
BOOST_AUTO_TEST_SUITE_END() // Examples
BOOST_AUTO_TEST_SUITE_END() // TestingTests
