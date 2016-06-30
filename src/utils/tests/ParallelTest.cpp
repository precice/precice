#include "ParallelTest.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"
#include "utils/Helpers.hpp"
#include <string>
#include <sstream>

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::utils::tests::ParallelTest)

namespace precice {
namespace utils {
namespace tests {

logging::Logger ParallelTest:: _log ( "precice::utils::tests::ParallelTest" );

ParallelTest:: ParallelTest ()
:
  TestCase ( "utils::tests::ParallelTest" )
{}

void ParallelTest:: run ()
{
# ifndef PRECICE_NO_MPI
  preciceTrace ( "run()" );
  typedef Parallel Par;
  if ( Par::getCommunicatorSize() >= 3 ){
    std::vector<int> ranksWanted;
    ranksWanted += 0, 1, 2;
    MPI_Comm comm = Par::getRestrictedCommunicator(ranksWanted);
    if ( Par::getProcessRank() <= 2 ){
      Par::setGlobalCommunicator(comm);
      std::string group;
      int rank = Par::getProcessRank();
      if ( (rank == 0) || (rank == 1) ){
        group = "GroupOne";
      }
      else {
        assertion ( rank == 2, rank );
        group = "GroupTwo";
      }
      Par::splitCommunicator( group );

      const std::vector<Par::AccessorGroup>& groups = Par::getAccessorGroups();
      validateEquals ( groups.size(), 2 );
      validateEquals ( groups[0].id, 0 );
      validateEquals ( groups[1].id, 1 );
      validateEquals ( groups[0].name, std::string("GroupOne") );
      validateEquals ( groups[1].name, std::string("GroupTwo") );
      validateEquals ( groups[0].leaderRank, 0 );
      validateEquals ( groups[1].leaderRank, 2 );
      validateEquals ( groups[0].size, 2 );
      validateEquals ( groups[1].size, 1 );

      Par::setGlobalCommunicator ( Par::getCommunicatorWorld() );
    }
  }
# endif // not PRECICE_NO_MPI
}


}}} // namespace precice, utils, tests
