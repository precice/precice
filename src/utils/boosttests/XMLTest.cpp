#include <string>
#include "testing/Testing.hpp"
#include "utils/xml/ValidatorEquals.hpp"
#include "utils/xml/ValidatorOr.hpp"
#include "utils/xml/XMLAttribute.hpp"
#include "utils/xml/XMLTag.hpp"

using namespace precice::utils;

BOOST_AUTO_TEST_SUITE(UtilsTests)
BOOST_AUTO_TEST_SUITE(XML)

struct CallbackHost : public XMLTag::Listener {
  Eigen::VectorXd eigenVectorXd;

  void xmlTagCallback(XMLTag &callingTag)
  {
    if (callingTag.getName() == "test-eigen-vectorxd-attributes") {
      eigenVectorXd = callingTag.getEigenVectorXdAttributeValue("value", 3);
    }
  }

  void xmlEndTagCallback(XMLTag &callingTag)
  {
    std::ignore = callingTag;
  }
};

BOOST_AUTO_TEST_CASE(AttributeConcatenation)
{
  std::string filename(Globals::getPathToSources() + "/utils/boosttests/config_xmltest_concatenation.xml");

  CallbackHost cb;
  XMLTag rootTag(cb, "configuration", XMLTag::OCCUR_ONCE);
  XMLTag testcaseTag(cb, "test-attribute-concatenation", XMLTag::OCCUR_ONCE);
  XMLTag testTag(cb, "test", XMLTag::OCCUR_ONCE_OR_MORE);

  XMLAttribute<std::string> attr("attribute");
  ValidatorEquals<std::string> equalsOne("value-one");
  ValidatorEquals<std::string> equalsTwo("value-two");

  ValidatorEquals<std::string> equalsThree("value-three");
  attr.setValidator(equalsOne || equalsTwo || equalsThree);
  testTag.addAttribute(attr);

  testcaseTag.addSubtag(testTag);
  rootTag.addSubtag(testcaseTag);

  configure(rootTag, filename);
}

BOOST_AUTO_TEST_CASE(VectorAttributes)
{
  std::string filename(Globals::getPathToSources() + "/utils/boosttests/config_xmltest_vectorattributes.xml");

  CallbackHost cb;
  XMLTag rootTag(cb, "configuration", XMLTag::OCCUR_ONCE);
  XMLTag testTagEigenXd(cb, "test-eigen-vectorxd-attributes", XMLTag::OCCUR_ONCE);

  XMLAttribute<Eigen::VectorXd> attrEigenXd("value");
  testTagEigenXd.addAttribute(attrEigenXd);
  rootTag.addSubtag(testTagEigenXd);

  configure(rootTag, filename);
  BOOST_TEST(cb.eigenVectorXd(0) == 3.0);
  BOOST_TEST(cb.eigenVectorXd(1) == 2.0);
  BOOST_TEST(cb.eigenVectorXd(2) == 1.0);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
