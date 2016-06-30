#include "XMLTest.hpp"
#include "utils/xml/XMLTag.hpp"
#include "utils/xml/XMLAttribute.hpp"
#include "utils/xml/ValidatorEquals.hpp"
#include "utils/xml/ValidatorOr.hpp"
#include "utils/Parallel.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Globals.hpp"
#include <string>

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::utils::tests::XMLTest)

namespace precice {
namespace utils {
namespace tests {

tarch::logging::Log XMLTest:: _log("precice::utils::XMLTest");

XMLTest:: XMLTest()
:
  tarch::tests::TestCase("precice::utils::tests::XMLTest"),
  _testDirectory(""),
  _vector2D(0.0),
  _vector3D(0.0),
  _dynVector()
{}

void XMLTest:: setUp()
{
  _testDirectory = Globals::getPathToSources() + "/utils/tests/";
}

//bool XMLTest:: initiateXMLFileReading
//(
//  XMLTag::XMLReader* xmlReader,
//  const std::string& filename )
//{
//  preciceTrace("initiateXMLFileReading()");
//  bool fileRead = true;
//  xmlReader = tarch::irr::io::createIrrXMLReader(filename.c_str());
//  if ((xmlReader==0)
//      || (not xmlReader->read())
//      || (xmlReader->getNodeType() == tarch::irr::io::EXN_NONE))
//  {
//     preciceError("initiateXMLFileReading()",
//                  "Was not able to read input file " + filename);
//     fileRead = false;
//  }
//  return fileRead;
//}

void XMLTest:: run()
{
  PRECICE_MASTER_ONLY {
    testMethod(testAttributeConcatenation);
    testMethod(testVectorAttributes);
    //testMethod(testNestedTags);
  }
}

void XMLTest:: testAttributeConcatenation()
{
  preciceTrace("testAttributeConcatenation()");
  std::string filename(_testDirectory + "config_xmltest_concatenation.xml");

  XMLTag rootTag(*this, "configuration", XMLTag::OCCUR_ONCE);
  XMLTag testcaseTag(*this, "test-attribute-concatenation", XMLTag::OCCUR_ONCE);
  XMLTag testTag(*this, "test", XMLTag::OCCUR_ONCE_OR_MORE);

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

void XMLTest:: testVectorAttributes()
{
  preciceTrace ( "testVectorAttributes()" );
  std::string filename (_testDirectory + "config_xmltest_vectorattributes.xml");

  XMLTag rootTag(*this, "configuration", XMLTag::OCCUR_ONCE);
  XMLTag testTag2D(*this, "test-vector2d-attributes", XMLTag::OCCUR_ONCE);
  XMLTag testTag3D(*this, "test-vector3d-attributes", XMLTag::OCCUR_ONCE);

  XMLAttribute<utils::Vector2D> attr2D("value");
  XMLAttribute<utils::Vector3D> attr3D("value");
  testTag2D.addAttribute(attr2D);
  testTag3D.addAttribute(attr3D);
  rootTag.addSubtag(testTag2D);
  rootTag.addSubtag(testTag3D);

  _vector2D = utils::Vector2D(0.0);
  _vector3D = utils::Vector3D(0.0);

  configure(rootTag, filename);
  validateEquals(_vector2D(0), 1.0);
  validateEquals(_vector2D(1), 2.0);
  validateEquals(_vector3D(0), 1.0);
  validateEquals(_vector3D(1), 2.0);
  validateEquals(_vector3D(2), 3.0);
}

void XMLTest:: xmlTagCallback
(
  XMLTag& callingTag )
{
  if (callingTag.getName() == "test-vector2d-attributes") {
    _vector2D = callingTag.getVector2DAttributeValue("value");
  }
  else if (callingTag.getName() == "test-vector3d-attributes") {
    _vector3D = callingTag.getVector3DAttributeValue("value");
  }
}

//void XMLTest:: testNestedTags()
//{
//  preciceTrace ( "testNestedTags()" );
//
//  std::string filename ( _testDirectory + "xmltest-nesting-config.xml" );
//  XMLTag tagNested(*this, "nested-tag", XMLTag::OCCUR_ARBITRARY_NESTED);
//  XMLTag subtag(*this, "subtag", XMLTag::OCCUR_ONCE);
//  XMLAttribute<int> attr("attribute");
//  subtag.addAttribute(attr);
//  tagNested.addSubtag(subtag);
//  bool success = configure(tagNested, filename);
//  validate(success);
//}

}}} // namespace precice, utils, tests
