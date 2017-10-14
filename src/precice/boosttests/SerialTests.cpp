#ifndef PRECICE_NO_MPI
#include "testing/Testing.hpp"

#include "precice/impl/SolverInterfaceImpl.hpp"
#include "precice/SolverInterface.hpp"
#include "precice/Constants.hpp"
#include "utils/Globals.hpp"
#include <fstream>

using namespace precice;

namespace precice {
extern bool testMode;
}


struct SerialTestFixture {
  SerialTestFixture()
  {
    precice::testMode = true;
  }

  ~SerialTestFixture()
  {
  }
};



BOOST_AUTO_TEST_SUITE(PreciceTests)
BOOST_FIXTURE_TEST_SUITE(Serial, SerialTestFixture)

// port all t-arch coupling mode tests here

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
#endif // PRECICE_NO_MPI
