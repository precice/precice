// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_TESTS_SERIALIMPLICITCOUPLINGSCHEMETEST_HPP_
#define PRECICE_CPLSCHEME_TESTS_SERIALIMPLICITCOUPLINGSCHEMETEST_HPP_

#include "com/Communication.hpp"
#include "m2n/M2N.hpp"
#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"
#include "utils/xml/XMLTag.hpp"
#include <string>
#include <vector>

namespace precice {
   namespace cplscheme {
      class CouplingScheme;
   }
   namespace mesh {
      class MeshConfiguration;
   }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace cplscheme {
namespace tests {

/**
 * @brief Tests class SerialImplicitCouplingScheme.
 */
class SerialImplicitCouplingSchemeTest : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor.
   */
  SerialImplicitCouplingSchemeTest();

  /**
   * @brief Destructor.
   */
  virtual ~SerialImplicitCouplingSchemeTest() {}

  /**
   * @brief Sets path to test directory.
   */
  virtual void setUp();

  /**
   * @brief Calls all test methods.
   */
  virtual void run();

private:

  // @brief Logging device.
  static tarch::logging::Log _log;

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
   * @brief Tests method extrapolateData.
   */
  void testExtrapolateData();

  /**
   * @brief Cpl. with absolute convergence measure and synchronized participants.
   *
   * ! General description
   * All objects are created and configured manually. A simple mesh consisting
   * of one vertex with two data is created. Each of the two participants is
   * writing one of the data and receiving the other data. No subcycling is
   * used by the participants.
   *
   * ! Test specific description
   * This test uses an absolute convergence measure for each data. Data values
   * are written with decreasing magnitude, which leads to convergence after
   * some coupling iterations. In this test the convergence behavior of the
   * participants is synchronized, i.e., they always perform the same amount
   * of iterations.
   */
  void testAbsConvergenceMeasureSynchronized();

  /**
   * @brief As testAbsConvergenceMeasureSynchronized, but configured from XML file.
   */
  void testConfiguredAbsConvergenceMeasureSynchronized();

  /**
   * @brief Minimum-iteration convergence measure, synchronized convergence.
   *
   * ! General description
   * The general description is equal to the one in
   * testAbsoluteConvergenceMeasure().
   *
   * ! Test specific description
   * This test uses a minimum convergence measure for each data. The
   * convergence behavior of the two participants is synchronized by giving
   * the same convergence measure for the data of each participant.
   */
  void testMinIterConvergenceMeasureSynchronized();

  /**
   * @brief Minimum-iteration convergence measure, assynchroneuous convergence.
   *
   * ! General description
   * The general description is equal to the one in
   * testAbsoluteConvergenceMeasure().
   *
   * ! Test specific description
   * This test uses a minimum convergence measure for each data. The
   * convergence behavior of the two participants is not synchronized, since
   * each participant uses a different number of minimum convergence cycles
   * to achieve convergence.
   */
  //   void testMinIterConvergenceMeasureAsync ();

  /**
   * @brief As testMinIterConvergenceMeasureAsync, but configured from XML file.
   */
  void testConfiguredMinIterConvergenceMeasureAsync();

  /**
   * @brief Performs tests described in testAbsoluteConvergenceMeasure().
   */
  void runCoupling (
    CouplingScheme&                cplScheme,
    const std::string&             nameParticipant,
    const mesh::MeshConfiguration& meshConfiguration,
    const std::vector<int>&        validIterations );

  /**
   * @brief Minimum-iteration convergence measure w. subcycling, synchronized.
   *
   * ! General description
   * All objects are created and configured manually. A simple mesh consisting
   * of one vertex with two data is created. Each of the two participants is
   * writing one of the data and receiving the other data. Subcycling with
   * different number of subcycles is applied by both participants.
   *
   * ! Test specific description
   * This test uses a minimum convergence measure for each data. The
   * convergence behavior of the two participants is synchronized by giving
   * the same convergence measure for the data of each participant.
   */
  void testMinIterConvergenceMeasureSynchronizedWithSubcycling();

  /**
   * @brief Tests the initialize data functionality.
   */
  void testInitializeData();

  /**
   * @brief Performs and validates implicit coupled simulation with subcycling.
   */
  void runCouplingWithSubcycling (
    CouplingScheme&                cplScheme,
    const std::string&             nameParticipant,
    const mesh::MeshConfiguration& meshConfig,
    const std::vector<int>&        validIterations );

  /**
   * @brief Performs and validates coupling with (first) participant setting dt.
   */
//  void runCouplingWithParticipantSetsDt (
//    CouplingScheme&                cplScheme,
//    const std::string&             nameParticipant,
//    const mesh::MeshConfiguration& meshConfig,
//    const std::vector<int>&        validIterations );

  void connect (
    const std::string&     participant0,
    const std::string&     participant1,
    const std::string&     localParticipant,
    m2n::M2N::SharedPointer& communication ) const;

# endif // not PRECICE_NO_MPI

};

}}} // namespace precice, cplscheme, tests

#endif /* PRECICE_CPLSCHEME_TESTS_SERIALIMPLICITCOUPLINGSCHEMETEST_HPP_ */
