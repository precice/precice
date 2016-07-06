#ifndef _TARCH_TESTS_TESTCASE_FACTORY_H_
#define _TARCH_TESTS_TESTCASE_FACTORY_H_

#include <string>

#include "tarch/tests/TestCase.h"


namespace tarch {
  namespace tests {
    template <class TestName>
    class TestCaseFactory;
  }
}


template <class TestName>
class tarch::tests::TestCaseFactory {
  private:
    TestCase* _test;
  public:
    enum Type {
      IntegrationTest,
      UnitTest
    };

    TestCaseFactory(Type type, const std::string& fullQualifiedTestName);

    ~TestCaseFactory();
};


#include "tarch/tests/TestCaseFactory.cpph"

#endif
