#ifndef PRECICE_NO_MPI

#include "testing/Testing.hpp"

#include <precice/precice.hpp>
#include <vector>

BOOST_AUTO_TEST_SUITE(Integration)
BOOST_AUTO_TEST_SUITE(Fundamental)
BOOST_AUTO_TEST_SUITE(JustInTime)
BOOST_AUTO_TEST_SUITE(Read)

PRECICE_TEST_SETUP("Receiver"_on(1_rank))
BOOST_AUTO_TEST_CASE(BeforeInitialize)
{
  PRECICE_TEST();

  // Set up Participant
  precice::Participant p(context.name, context.config(), 0, 1);

  std::vector<double> readAt = {0.25, 0.0, 0.75, 0.0};
  std::vector<double> readData(2, 0.0);
  BOOST_CHECK_EXCEPTION(p.mapAndReadData("M", "D", readAt, 0.0, readData), precice::Error, ::precice::testing::errorContains("before initialize"));
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

#endif // PRECICE_NO_MPI
