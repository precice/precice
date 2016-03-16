// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_QUERY_FINDCLOSESTVERTEXVISITORTEST_HPP_
#define PRECICE_QUERY_FINDCLOSESTVERTEXVISITORTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace query {
namespace tests {


class FindClosestVertexVisitorTest : public tarch::tests::TestCase
{
public:

   FindClosestVertexVisitorTest();

   virtual ~FindClosestVertexVisitorTest() {}

   virtual void setUp() {}

   virtual void run();

private:

   static logging::Logger _log;
};

}}} // namespace precice, query, tests

#endif /* PRECICE_QUERY_FINDCLOSESTVERTEXVISITORTEST_HPP_ */
