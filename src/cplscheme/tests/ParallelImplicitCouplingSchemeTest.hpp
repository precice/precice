#pragma once

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "m2n/SharedPointer.hpp"
#include <string>

namespace precice {
   namespace cplscheme {
      class CouplingScheme;
   }
   namespace mesh {
      class MeshConfiguration;
   }
}

namespace precice {
namespace cplscheme {
namespace tests {

/**
 * @brief Tests class ParallelImplicitCouplingSchemeTest.
 */
class ParallelImplicitCouplingSchemeTest : public tarch::tests::TestCase
{
public:

  ParallelImplicitCouplingSchemeTest();

  virtual ~ParallelImplicitCouplingSchemeTest() {}

  /**
   * @brief Sets path to test directory.
   */
  virtual void setUp();

  /**
   * @brief Calls all test methods.
   */
  virtual void run();

  typedef std::map<int,PtrCouplingData> DataMap;
  
private:

  // @brief Logging device.
  static logging::Logger _log;

  // @brief Holds file path to precice src directory.
  std::string _pathToTests;

  const std::string MY_WRITE_CHECKPOINT;

  const std::string MY_READ_CHECKPOINT;

  //utils::XMLTag _root;

# ifndef PRECICE_NO_MPI

  /**
   * @brief Tests reading of XML config with relaxed coupling iterations.
   */
  void testParseConfigurationWithRelaxation ();

  /**
   * @brief Tests the initialize data functionality.
   */
  void testInitializeData();

  /**
   * @brief Tests the correct postprocessing for VIQN-like vector data
   */
  void testVIQNPP();
  
  /**
   * @brief Tests the correct postprocessing for MVQN-like vector data
   */
  void testMVQNPP();

  void connect (
      const std::string&     participant0,
      const std::string&     participant1,
      const std::string&     localParticipant,
      m2n::PtrM2N&           communication ) const;

# endif // not PRECICE_NO_MPI
};

}}} // namespace precice, cplscheme, tests
