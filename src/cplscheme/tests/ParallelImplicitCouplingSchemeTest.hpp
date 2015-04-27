// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_TESTS_PARALLELIMPLICITCOUPLINGSCHEMETEST_HPP_
#define PRECICE_CPLSCHEME_TESTS_PARALLELIMPLICITCOUPLINGSCHEMETEST_HPP_

#include "com/Communication.hpp"
#include "m2n/M2N.hpp"
#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"
#include "utils/xml/XMLTag.hpp"
#include "cplscheme/SharedPointer.hpp"
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
 * @brief Tests class ParallelImplicitCouplingSchemeTest.
 */
class ParallelImplicitCouplingSchemeTest : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor.
   */
  ParallelImplicitCouplingSchemeTest();

  /**
   * @brief Destructor.
   */
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
  typedef tarch::la::DynamicColumnMatrix<double> DataMatrix;

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
      m2n::M2N::SharedPointer&           communication ) const;

# endif // not PRECICE_NO_MPI
};

}}} // namespace precice, cplscheme, tests

#endif /* PRECICE_CPLSCHEME_TESTS_PARALLELIMPLICITCOUPLINGSCHEMETEST_HPP_ */
