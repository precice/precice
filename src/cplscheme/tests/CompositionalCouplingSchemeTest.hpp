#pragma once

#include "cplscheme/SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"
#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"
#include "m2n/SharedPointer.hpp"
#include <string>

namespace precice {
  namespace mesh {
    class MeshConfiguration;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace cplscheme {
namespace tests {

/// Provides unit tests for class CompositionalCouplingScheme.
class CompositionalCouplingSchemeTest : public tarch::tests::TestCase
{
public:

  CompositionalCouplingSchemeTest();

  virtual ~CompositionalCouplingSchemeTest() {}

  /// Sets path to test directory.
  virtual void setUp();

  /**
   * @brief Runs all tests.
   */
  virtual void run();

private:

  // @brief Logging device.
  static logging::Logger _log;

  // @brief Path to src directory.
  std::string _pathToTests;

  /**
   * @brief Runs different explicit/implicit compositions with dummy.
   */
  void testDummySchemeCompositions();

# ifndef PRECICE_NO_MPI

  /**
   * @brief Runs a three solver composition of explicit schemes.
   */
  void testExplicitSchemeComposition1();

  /**
   * @brief Runs a three solver composition of implicit schemes.
   */
  void testImplicitSchemeComposition();

  /**
   * @brief Runs a three solver composition, S1 <-impl.-> S2 <-expl.-> S3.
   */
  void testImplicitExplicitSchemeComposition();

  /**
   * @brief Runs a three solver composition, S1 <-expl.-> S2 <-impl.-> S3.
   */
  void testExplicitImplicitSchemeComposition();

  /**
   * @brief Setup three solver coupling using XML-configuration.
   */
  void setupAndRunThreeSolverCoupling(const std::string& configFilename);

  /**
   * @brief Takes a configured coupling scheme and performs explicit coupling.
   *
   * @param[in] participantName Either "participant0" or "participant1".
   */
  void runThreeSolverCoupling (
    PtrCouplingScheme          cplScheme0,
    const std::string&         participantName,
    mesh::PtrMeshConfiguration meshConfig );

  void connect (
    const std::string&     participant0,
    const std::string&     participant1,
    const std::string&     localParticipant,
    m2n::PtrM2N& communication ) const;


# endif // not PRECICE_NO_MPI
};

}}} // namespace precice, cplscheme, tests

