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

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(MissingRequiredAttribute)
{
  PRECICE_TEST();
  // A required attribute with no default value is absent from the XML file.
  // The error message should include the owning tag name so the user can
  // locate the problematic element without grepping their config manually.
  std::string filename(getPathToSources() + "/xml/tests/xmltest-missing-attr.xml");

  CallbackHost cb;
  XMLTag       rootTag(cb, "configuration", XMLTag::OCCUR_ONCE);
  XMLTag       testTag(cb, "test-tag", XMLTag::OCCUR_ONCE);

  // "name" is required (no default given)
  XMLAttribute<std::string> attr("name");
  testTag.addAttribute(attr);
  rootTag.addSubtag(testTag);

  BOOST_CHECK_EXCEPTION(
      configure(rootTag, ConfigurationContext{}, filename),
      ::precice::Error,
      ::precice::testing::errorContains("test-tag"));
}

PRECICE_TEST_SETUP(1_rank)
BOOST_AUTO_TEST_CASE(InvalidAttributeOption)
{
  PRECICE_TEST();
  // An attribute is present but its value is not in the declared option set.
  // The error message should include the owning tag name so the user knows
  // which element to fix.
  std::string filename(getPathToSources() + "/xml/tests/xmltest-invalid-attr-option.xml");

  CallbackHost cb;
  XMLTag       rootTag(cb, "configuration", XMLTag::OCCUR_ONCE);
  XMLTag       testTag(cb, "test-tag", XMLTag::OCCUR_ONCE);

  // "mode" only accepts "fast" or "slow"
  auto attr = makeXMLAttribute("mode", "fast").setOptions({"fast", "slow"});
  testTag.addAttribute(attr);
  rootTag.addSubtag(testTag);

  BOOST_CHECK_EXCEPTION(
      configure(rootTag, ConfigurationContext{}, filename),
      ::precice::Error,
      ::precice::testing::errorContains("test-tag"));
}

BOOST_AUTO_TEST_SUITE_END()
