#include <Eigen/Core>
#include <optional>
#include <string>
#include <tuple>
#include <vector>
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
      eigenVectorXd = callingTag.getEigenVectorXdAttributeValue("value");
    }
  }

  void xmlEndTagCallback(const ConfigurationContext &context, XMLTag &callingTag) override
  {
    std::ignore = callingTag;
  }
};

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(AttributeConcatenation)
{
  PRECICE_TEST();
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

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(VectorAttributes)
{
  PRECICE_TEST();
  std::string filename(getPathToSources() + "/xml/tests/config_xmltest_vectorattributes.xml");

  CallbackHost cb;
  XMLTag       rootTag(cb, "configuration", XMLTag::OCCUR_ONCE);
  XMLTag       testTagEigenXd(cb, "test-eigen-vectorxd-attributes", XMLTag::OCCUR_ONCE);

  XMLAttribute<Eigen::VectorXd> attrEigenXd("value");
  testTagEigenXd.addAttribute(attrEigenXd);
  rootTag.addSubtag(testTagEigenXd);

  configure(rootTag, ConfigurationContext{}, filename);
  BOOST_TEST(cb.eigenVectorXd.size() == 3);
  BOOST_TEST(cb.eigenVectorXd(0) == 3.0);
  BOOST_TEST(cb.eigenVectorXd(1) == 2.0);
  BOOST_TEST(cb.eigenVectorXd(2) == 1.0);
}

struct OptionalAttributeHost : public XMLTag::Listener {
  std::vector<std::optional<int>> values;

  void xmlTagCallback(const ConfigurationContext &, XMLTag &callingTag) override
  {
    if (callingTag.getName() == "test-optional-attribute") {
      values.push_back(callingTag.getOptionalIntAttributeValue("optional-value"));
    }
  }

  void xmlEndTagCallback(const ConfigurationContext &, XMLTag &) override {}
};

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(OptionalAttributeWithoutDefault)
{
  PRECICE_TEST();

  std::string filename(getPathToSources() + "/xml/tests/config_xmltest_optional_attribute.xml");

  OptionalAttributeHost cb;
  XMLTag                rootTag(cb, "configuration", XMLTag::OCCUR_ONCE);
  XMLTag                testTag(cb, "test-optional-attribute", XMLTag::OCCUR_ONCE_OR_MORE);

  XMLAttribute<int> attrOptional("optional-value");
  attrOptional.setOptional();
  testTag.addAttribute(attrOptional);

  rootTag.addSubtag(testTag);

  configure(rootTag, ConfigurationContext{}, filename);

  BOOST_TEST(cb.values.size() == 2u);
  BOOST_TEST(cb.values[0].has_value());
  BOOST_TEST(cb.values[0].value() == 42);
  BOOST_TEST(!cb.values[1].has_value());
}

BOOST_AUTO_TEST_SUITE_END()
