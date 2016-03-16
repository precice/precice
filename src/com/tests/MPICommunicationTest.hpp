// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_COM_TESTS_COMMUNICATIONMPITEST_HPP_
#define PRECICE_COM_TESTS_COMMUNICATIONMPITEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace com {
namespace tests {

/**
 * @brief Provides tests for class MPIPortsCommunication.
 *
 * Since MPICommunication is abstract, MPIPortsCommunication is used wihtout
 * establishing a porst connection. Thus, only local communication is tested.
 */
class MPICommunicationTest : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor.
   */
   MPICommunicationTest();

   /**
    * @brief Destructor, empty.
    */
   virtual ~MPICommunicationTest() {}

   /**
    * @brief Empty.
    */
   virtual void setUp() {}

   /**
    * @brief Runs all tests.
    */
   virtual void run();

private:

   // @brief Logging device.
   static logging::Logger _log;

#  ifndef PRECICE_NO_MPI

   void testSendAndReceiveString();

   void testSendAndReceiveVector();

   void testSendAndReceiveInteger();

#  endif // not PRECICE_NO_MPI
};

}}} // namespace precice, com, tests

#endif /* PRECICE_COM_TESTS_COMMUNICATIONMPITEST_HPP_ */
