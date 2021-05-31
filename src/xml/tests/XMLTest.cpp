#include <Eigen/Core>
#include <string>
#include <tuple>
#include "math/constants.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "xml/ValueParser.hpp"
#include "xml/XMLAttribute.hpp"
#include "xml/XMLTag.hpp"

using namespace precice;
using namespace precice::xml;
using precice::testing::getPathToSources;

BOOST_AUTO_TEST_SUITE(XML)

struct CallbackHost : public XMLTag::Listener {
  Eigen::VectorXd eigenVectorXd;

  void xmlTagCallback(const ConfigurationContext &context, XMLTag &callingTag) override
  {
    if (callingTag.getName() == "test-eigen-vectorxd-attributes") {
      eigenVectorXd = callingTag.getEigenVectorXdAttributeValue("value", 3);
    }
  }

  void xmlEndTagCallback(const ConfigurationContext &context, XMLTag &callingTag) override
  {
    std::ignore = callingTag;
  }
};

BOOST_AUTO_TEST_CASE(AttributeConcatenation)
{
  PRECICE_TEST(1_rank);
  std::string filename(getPathToSources() + "/xml/tests/config_xmltest_concatenation.xml");

  CallbackHost cb;
  XMLTag       rootTag(cb, "configuration", XMLTag::OCCUR_ONCE);
  XMLTag       testcaseTag(cb, "test-attribute-concatenation", XMLTag::OCCUR_ONCE);
  XMLTag       testTag(cb, "test", XMLTag::OCCUR_ONCE_OR_MORE);

  auto attr = makeXMLAttribute("attribute", "").setOptions({"value-one", "value-two", "value-three"});
  testTag.addAttribute(attr);

  testcaseTag.addSubtag(testTag);
  rootTag.addSubtag(testcaseTag);

  configure(rootTag, ConfigurationContext{}, filename);
}

BOOST_AUTO_TEST_CASE(VectorAttributes)
{
  PRECICE_TEST(1_rank);
  std::string filename(getPathToSources() + "/xml/tests/config_xmltest_vectorattributes.xml");

  CallbackHost cb;
  XMLTag       rootTag(cb, "configuration", XMLTag::OCCUR_ONCE);
  XMLTag       testTagEigenXd(cb, "test-eigen-vectorxd-attributes", XMLTag::OCCUR_ONCE);

  XMLAttribute<Eigen::VectorXd> attrEigenXd("value");
  testTagEigenXd.addAttribute(attrEigenXd);
  rootTag.addSubtag(testTagEigenXd);

  configure(rootTag, ConfigurationContext{}, filename);
  BOOST_TEST(cb.eigenVectorXd(0) == 3.0);
  BOOST_TEST(cb.eigenVectorXd(1) == 2.0);
  BOOST_TEST(cb.eigenVectorXd(2) == 1.0);
}

BOOST_AUTO_TEST_SUITE_END()
