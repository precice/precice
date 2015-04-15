// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_TESTS_COMPOSITIONALCOUPLINGSCHEMETEST_HPP_
#define PRECICE_CPLSCHEME_TESTS_COMPOSITIONALCOUPLINGSCHEMETEST_HPP_

#include "cplscheme/SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"
#include "com/Communication.hpp"
#include "m2n/M2N.hpp"
#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"
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

/**
 * @brief Provides unit tests for class CompositionalCouplingScheme.
 */
class CompositionalCouplingSchemeTest : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor.
   */
  CompositionalCouplingSchemeTest();

  /**
   * @brief Destructor, empty.
   */
  virtual ~CompositionalCouplingSchemeTest() {}

  /**
   * @brief Sets path to test directory.
   */
  virtual void setUp();

  /**
   * @brief Runs all tests.
   */
  virtual void run();

private:

  // @brief Logging device.
  static tarch::logging::Log _log;

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
   * @param participantName [IN] Either "participant0" or "participant1".
   */
  void runThreeSolverCoupling (
    PtrCouplingScheme          cplScheme0,
    const std::string&         participantName,
    mesh::PtrMeshConfiguration meshConfig );

  void connect (
    const std::string&     participant0,
    const std::string&     participant1,
    const std::string&     localParticipant,
    m2n::M2N::SharedPointer& communication ) const;


# endif // not PRECICE_NO_MPI
};

}}} // namespace precice, cplscheme, tests

#endif /* PRECICE_CPLSCHEME_TESTS_COMPOSITIONALCOUPLINGSCHEMETEST_HPP_ */
