#include "ManageUniqueIDsTest.hpp"
#include "utils/ManageUniqueIDs.hpp"
#include "utils/Parallel.hpp"
#include "utils/Helpers.hpp"
#include "utils/Globals.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::utils::tests::ManageUniqueIDsTest)

namespace precice {
namespace utils {
namespace tests {

tarch::logging::Log ManageUniqueIDsTest::
   _log ( "precice::utils::tests::ManageUniqueIDsTest" );

ManageUniqueIDsTest:: ManageUniqueIDsTest ()
:
   TestCase ("utils::ManageUniqueIDsTest")
{}

void ManageUniqueIDsTest:: run ()
{
   PRECICE_MASTER_ONLY {
      testMethod(testUniqueIDs);
   }
}

void ManageUniqueIDsTest:: testUniqueIDs ()
{
  preciceTrace ( "testUniqueIDs()" );
  ManageUniqueIDs uniqueIDs;
  int id = uniqueIDs.getFreeID ();
  validateEquals (id, 0);
  id = uniqueIDs.getFreeID ();
  validateEquals (id, 1);
  bool success = uniqueIDs.insertID (2);
  validate (success);
  id = uniqueIDs.getFreeID ();
  validateEquals (id, 3);
}

}}} // namespace precice, utils, tests
