#include <Eigen/Core>
#include <string>
#include <tuple>
#include "math/constants.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "xml/ConfigParser.hpp"
#include "xml/ValueParser.hpp"
#include "xml/XMLAttribute.hpp"
#include "xml/XMLTag.hpp"

using namespace precice;
using namespace precice::xml;
using precice::testing::getPathToSources;

BOOST_AUTO_TEST_SUITE(XML)

struct CallbackHostAttr : public XMLTag::Listener {
  Eigen::VectorXd eigenValue;
  double          doubleValue;
  int             intValue;
  bool            boolValue;
  std::string     stringValue;

  void xmlTagCallback(const ConfigurationContext &context, XMLTag &callingTag) override
  {
    if (callingTag.getName() == "test-double") {
      doubleValue = callingTag.getDoubleAttributeValue("attribute");
    }

    if (callingTag.getName() == "test-eigen") {
      eigenValue = callingTag.getEigenVectorXdAttributeValue("attribute", 3);
    }

    if (callingTag.getName() == "test-int") {
      intValue = callingTag.getIntAttributeValue("attribute");
    }

    if (callingTag.getName() == "test-bool") {
      boolValue = callingTag.getBooleanAttributeValue("attribute");
    }

    if (callingTag.getName() == "test-string") {
      stringValue = callingTag.getStringAttributeValue("text");
    }
  }

  void xmlEndTagCallback(const ConfigurationContext &context, XMLTag &callingTag) override
  {
    std::ignore = callingTag;
  }
};

BOOST_AUTO_TEST_CASE(AttributeTypeTest)
{
  PRECICE_TEST(1_rank);
  std::string filename(getPathToSources() + "/xml/tests/xmlparser_test.xml");

  CallbackHostAttr cb;
  XMLTag           rootTag(cb, "configuration", XMLTag::OCCUR_ONCE);
  XMLTag           testcaseTag(cb, "test-config", XMLTag::OCCUR_ONCE);

  XMLTag doubleTag(cb, "test-double", XMLTag::OCCUR_ONCE_OR_MORE);
  XMLTag eigenTag(cb, "test-eigen", XMLTag::OCCUR_ONCE_OR_MORE);
  XMLTag intTag(cb, "test-int", XMLTag::OCCUR_ONCE_OR_MORE);
  XMLTag stringTag(cb, "test-string", XMLTag::OCCUR_ONCE_OR_MORE);
  XMLTag boolTag(cb, "test-bool", XMLTag::OCCUR_ONCE_OR_MORE);

  XMLAttribute<double>          doubleAttr("attribute");
  XMLAttribute<Eigen::VectorXd> eigenAttr("attribute");
  XMLAttribute<int>             intAttr("attribute");
  XMLAttribute<bool>            boolAttr("attribute");
  XMLAttribute<std::string>     stringAttr("text");

  doubleTag.addAttribute(doubleAttr);
  eigenTag.addAttribute(eigenAttr);
  intTag.addAttribute(intAttr);
  boolTag.addAttribute(boolAttr);
  stringTag.addAttribute(stringAttr);

  testcaseTag.addSubtag(doubleTag);
  testcaseTag.addSubtag(eigenTag);
  testcaseTag.addSubtag(intTag);
  testcaseTag.addSubtag(boolTag);
  testcaseTag.addSubtag(stringTag);

  rootTag.addSubtag(testcaseTag);

  configure(rootTag, ConfigurationContext{}, filename);

  BOOST_TEST(cb.boolValue == true);
  BOOST_TEST(cb.doubleValue == 3.1);
  BOOST_TEST(cb.intValue == 4);
  BOOST_TEST(cb.stringValue == "Hello World");

  BOOST_TEST(cb.eigenValue(0) == 3.0);
  BOOST_TEST(cb.eigenValue(1) == 2.0);
  BOOST_TEST(cb.eigenValue(2) == 1.0);
}

BOOST_AUTO_TEST_CASE(OccurenceTest)
{
  PRECICE_TEST(1_rank);
  std::string filename(getPathToSources() + "/xml/tests/xmlparser_occtest.xml");

  CallbackHostAttr cb;
  XMLTag           rootTag(cb, "configuration", XMLTag::OCCUR_ONCE);
  XMLTag           testcaseTag(cb, "test-config", XMLTag::OCCUR_ONCE);

  XMLTag occ2(cb, "test-occ_once_or_more-once", XMLTag::OCCUR_ONCE_OR_MORE);
  XMLTag occ3(cb, "test-occ_arbitrary", XMLTag::OCCUR_ARBITRARY);
  XMLTag occ4(cb, "test-occ_not_or_once", XMLTag::OCCUR_NOT_OR_ONCE);

  XMLTag occ5(cb, "test-occ_once_or_more-more", XMLTag::OCCUR_ONCE_OR_MORE);
  XMLTag occ6(cb, "test-occ_arbitrary-opt", XMLTag::OCCUR_ARBITRARY);
  XMLTag occ7(cb, "test-occ_not_or_once-opt", XMLTag::OCCUR_NOT_OR_ONCE);

  testcaseTag.addSubtag(occ2);
  testcaseTag.addSubtag(occ3);
  testcaseTag.addSubtag(occ4);

  testcaseTag.addSubtag(occ5);
  testcaseTag.addSubtag(occ6);
  testcaseTag.addSubtag(occ7);

  rootTag.addSubtag(testcaseTag);

  configure(rootTag, ConfigurationContext{}, filename);
}

BOOST_AUTO_TEST_CASE(NamespaceTest)
{
  PRECICE_TEST(1_rank);
  std::string filename(getPathToSources() + "/xml/tests/xmlparser_nstest.xml");

  CallbackHostAttr cb;
  XMLTag           rootTag(cb, "configuration", XMLTag::OCCUR_ONCE);
  XMLTag           testcaseTag(cb, "test-config", XMLTag::OCCUR_ONCE);

  XMLTag intTag(cb, "test-no_ns", XMLTag::OCCUR_ONCE_OR_MORE);
  XMLTag stringTag(cb, "test-ns", XMLTag::OCCUR_ONCE_OR_MORE, "ns");

  testcaseTag.addNamespace("ns");

  testcaseTag.addSubtag(intTag);
  testcaseTag.addSubtag(stringTag);

  rootTag.addSubtag(testcaseTag);

  configure(rootTag, ConfigurationContext{}, filename);
}

struct ContextListener : public XMLTag::Listener {
  ConfigurationContext startContext;
  ConfigurationContext endContext;

  void xmlTagCallback(const ConfigurationContext &context, XMLTag &callingTag)
  {
    startContext = context;
  }

  void xmlEndTagCallback(const ConfigurationContext &context, XMLTag &callingTag)
  {
    endContext = context;
  }
};

BOOST_AUTO_TEST_CASE(Context)
{
  PRECICE_TEST(1_rank);
  std::string filename(getPathToSources() + "/xml/tests/config_xmltest_context.xml");

  ContextListener      cl;
  XMLTag               rootTag(cl, "configuration", XMLTag::OCCUR_ONCE);
  ConfigurationContext ccontext{"test", 12, 32};
  configure(rootTag, ccontext, filename);
  BOOST_TEST(cl.startContext.name == "test");
  BOOST_TEST(cl.startContext.rank == 12);
  BOOST_TEST(cl.startContext.size == 32);
  BOOST_TEST(cl.endContext.name == "test");
  BOOST_TEST(cl.endContext.rank == 12);
  BOOST_TEST(cl.endContext.size == 32);
}

BOOST_AUTO_TEST_CASE(Decode)
{
  PRECICE_TEST(1_rank);

  BOOST_TEST(decodeXML("Less than &lt; test") == "Less than < test");
  BOOST_TEST(decodeXML("Greater than &gt; test") == "Greater than > test");
  BOOST_TEST(decodeXML("Ampersand &amp; test") == "Ampersand & test");
  BOOST_TEST(decodeXML("Quotation &quot; test") == "Quotation \" test");
  BOOST_TEST(decodeXML("Apostrophe &apos; test") == "Apostrophe ' test");

  BOOST_TEST(decodeXML("&quot; &lt; &gt; &gt; &lt; &amp; &quot; &amp; &apos;") == "\" < > > < & \" & '");
}

BOOST_AUTO_TEST_SUITE_END()
