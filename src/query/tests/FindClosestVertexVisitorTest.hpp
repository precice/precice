#ifndef PRECICE_QUERY_FINDCLOSESTVERTEXVISITORTEST_HPP_
#define PRECICE_QUERY_FINDCLOSESTVERTEXVISITORTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "tarch/logging/Log.h"

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

   static tarch::logging::Log _log;
};

}}} // namespace precice, query, tests

#endif /* PRECICE_QUERY_FINDCLOSESTVERTEXVISITORTEST_HPP_ */
