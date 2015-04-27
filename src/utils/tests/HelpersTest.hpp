// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_UTILS_TESTS_HELPERSTEST_HPP_
#define PRECICE_UTILS_TESTS_HELPERSTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

namespace precice {
namespace utils {
namespace tests {

/**
 * @brief Tests methods from file Helpers.
 */
class HelpersTest : public tarch::tests::TestCase
{
public:

   HelpersTest ();

   virtual void setUp () {}

   virtual void run ();

private:

    static tarch::logging::Log _log;

//    void testGetZero ();

    void testAppendTo ();

    void testOperatorPlusForVectors ();
};

}}} // namespace precice, utils, tests

#endif /* PRECICE_UTILS_TESTS_HELPERSTEST_HPP_ */
