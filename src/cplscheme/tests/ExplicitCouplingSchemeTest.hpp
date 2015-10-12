// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_TESTS_EXPLICITCOUPLINGSCHEMETEST_HPP_
#define PRECICE_CPLSCHEME_TESTS_EXPLICITCOUPLINGSCHEMETEST_HPP_

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
  namespace com {
    class CommunicationConfiguration;
  }
  namespace cplscheme {
    class CouplingScheme;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace cplscheme {
namespace tests {

class ExplicitCouplingSchemeTest : public tarch::tests::TestCase
{
public:

   /**
    * @brief Constructor.
    */
   ExplicitCouplingSchemeTest();

   /**
    * @brief Destructor, empty.
    */
   virtual ~ExplicitCouplingSchemeTest() {}

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

#  ifndef PRECICE_NO_MPI

   /**
    * @brief Runs a simple explicit coupling scheme with two participants.
    *
    * As mesh, only one vertex is created, that holds two data. Each participant
    * does send one of the data and receive the other one. No subcycling is
    * applied. The data values are kown apriori, in order to validate them.
    * The data values are set by the solvers and increased in every timestep,
    * in order to see, if a proper transmission of values occurs at every time-
    * step.
    */
   void testSimpleExplicitCoupling();

   /**
    * @brief Same as testSimpleExplicitCoupling, but configured from XML file.
    */
   void testConfiguredSimpleExplicitCoupling();

   /**
    * @brief Configured test with first participant setting timestep length.
    */
   void testExplicitCouplingFirstParticipantSetsDt();

   /**
    * @brief Test from XML configuration to test data initialization for serial
    * coupling.
    */
   void testSerialDataInitialization();

   /**
    * @brief Test from XML configuration to test data initialization for parallel
    * coupling.
    */
   void testParallelDataInitialization();

   /**
    * @brief Configured test with second participant setting timestep length.
    */
   //void testExplicitCouplingSecondParticipantSetsDt();

   /**
    * @brief Takes a configured coupling scheme and performs explicit coupling.
    *
    * @param cplScheme [IN/OUT] Reference to explicit coupling scheme
    * @param participantName [IN] Either "participant0" or "participant1".
    */
   void runSimpleExplicitCoupling (
      CouplingScheme &        cplScheme,
      const std::string  &            participantName,
      const mesh::MeshConfiguration & meshConfig );

   /**
    * @brief Runs an explicit coupling scheme with two participants, subcycling.
    *
    * As mesh, only one vertex is created, that holds two data. Each participant
    * does send one of the data and receive the other one. One of the
    * participants employs subcycling. The data values are kown apriori, in
    * order to validate them. The data values are set by the solvers and
    * increased in every timestep, in order to see, if a proper transmission of
    * values occurs at every timestep.
    */
   void testExplicitCouplingWithSubcycling ();

   /**
    * @brief Same as testExplicitCouplingWithSubcycling, but configured from XML file.
    */
   void testConfiguredExplicitCouplingWithSubcycling ();

   void runExplicitCouplingWithSubcycling (
      CouplingScheme &        cplScheme,
      const std::string &             participantName,
      const mesh::MeshConfiguration & meshConfig );

   void connect (
     const std::string &      participant0,
     const std::string &      participant1,
     const std::string &      localParticipant,
     m2n::M2N::SharedPointer & communication ) const;


#  endif // not PRECICE_NO_MPI
};

}}} // namespace precice, cplscheme, tests

#endif /* PRECICE_CPLSCHEME_TESTS_EXPLICITCOUPLINGSCHEMETEST_HPP_ */
