#include <algorithm>
#include <boost/test/framework.hpp>
#include <boost/test/tree/test_unit.hpp>
#include <precice/Exceptions.hpp>

#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <limits>
#include <regex>
#include <sstream>
#include <string>

#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xpath.h>

#include "logging/LogMacros.hpp"
#include "logging/Logger.hpp"
#include "testing/SourceLocation.hpp"
#include "testing/Testing.hpp"
#include "utils/ManageUniqueIDs.hpp"
#include "utils/assertion.hpp"

namespace precice::testing {

DataID operator""_dataID(unsigned long long n)
{
  PRECICE_ASSERT(n < std::numeric_limits<DataID>::max(), "DataID is too big");
  return static_cast<DataID>(n);
}

std::string getPathToRepository()
{
  precice::logging::Logger _log("testing");
  return precice::testing::SourceLocation;
}

std::string getPathToSources()
{
  return getPathToRepository() + "/src";
}

std::string getPathToTests()
{
  return getPathToRepository() + "/tests";
}

std::string getTestPath()
{
  const auto &cspan = boost::unit_test::framework::current_test_case().p_file_name;
  return {cspan.begin(), cspan.end()};
}

std::string getTestName()
{
  std::string name = boost::unit_test::framework::current_test_case().p_name;
  // check for data tests
  if (name[0] != '_') {
    return name;
  }
  auto parent = boost::unit_test::framework::current_test_case().p_parent_id;
  return boost::unit_test::framework::get<boost::unit_test::test_suite>(parent).p_name;
}

std::string getFullTestName()
{
  return boost::unit_test::framework::current_test_case().full_name();
}

TestSetup getTestSetup()
{
  if (auto setup = getTestSetupFor(boost::unit_test::framework::current_test_case());
      setup) {
    return *setup;
  }
  throw std::runtime_error{"Use PRECICE_TEST_SETUP(...) to define the setup directly before the BOOST_[AUTO_]TEST_CASE."};
}

std::optional<TestSetup> getTestSetupFor(const boost::unit_test::test_unit &tu)
{
  const auto &list = tu.p_decorators.get();
  for (const auto &iter : list) {
    auto ptr = dynamic_cast<precice::testing::precice_testsetup_fixture *>(iter.get());
    if (ptr) {
      return ptr->testSetup;
    }
  }
  return std::nullopt;
}

int nextMeshID()
{
  static utils::ManageUniqueIDs manager;
  return manager.getFreeID();
}

/// equals to be used in tests. Compares two std::vectors using a given tolerance. Prints both operands on failure
boost::test_tools::predicate_result equals(const std::vector<float> &VectorA,
                                           const std::vector<float> &VectorB,
                                           float                     tolerance)
{
  PRECICE_ASSERT(VectorA.size() == VectorB.size());
  Eigen::MatrixXd MatrixA(VectorA.size(), 1);
  std::copy(VectorA.begin(), VectorA.end(), MatrixA.data());
  Eigen::MatrixXd MatrixB(VectorB.size(), 1);
  std::copy(VectorB.begin(), VectorB.end(), MatrixB.data());
  return equals(MatrixA, MatrixB, tolerance);
}

boost::test_tools::predicate_result equals(const std::vector<double> &VectorA,
                                           const std::vector<double> &VectorB,
                                           double                     tolerance)
{
  PRECICE_ASSERT(VectorA.size() == VectorB.size());
  Eigen::MatrixXd MatrixA(VectorA.size(), 1);
  std::copy(VectorA.begin(), VectorA.end(), MatrixA.data());
  Eigen::MatrixXd MatrixB(VectorB.size(), 1);
  std::copy(VectorB.begin(), VectorB.end(), MatrixB.data());
  return equals(MatrixA, MatrixB, tolerance);
}

boost::test_tools::predicate_result equals(float a, float b, float tolerance)
{
  if (not math::equals(a, b, tolerance)) {
    boost::test_tools::predicate_result res(false);
    res.message() << "Not equal: " << a << "!=" << b;
    return res;
  }
  return true;
}

boost::test_tools::predicate_result equals(double a, double b, double tolerance)
{
  if (not math::equals(a, b, tolerance)) {
    boost::test_tools::predicate_result res(false);
    res.message() << "Not equal: " << a << "!=" << b;
    return res;
  }
  return true;
}

void expectFile(std::string_view name)
{
  BOOST_TEST(std::filesystem::is_regular_file(name), "File " << name << " is not a regular file or doesn't exist.");
}

void removeFile(std::string_view name)
{
  if (std::filesystem::exists(name)) {
    std::filesystem::remove(name);
  }
}

std::string readFile(std::string_view name)
{
  std::ifstream t{std::string(name)};
  std::stringstream buffer;
  buffer << t.rdbuf();
  return buffer.str();
}

static void xmlErrorSuppress(void*, const char*, ...) {}

bool validateXMLWellFormed(std::string_view filename)
{
  xmlSetGenericErrorFunc(nullptr, xmlErrorSuppress);
  xmlDocPtr doc = xmlReadFile(std::string(filename).c_str(), nullptr, 
                               XML_PARSE_NOWARNING | XML_PARSE_NOERROR);
  bool valid = (doc != nullptr);
  if (doc) xmlFreeDoc(doc);
  return valid;
}

static bool validateVTKFile(std::string_view filename, const char* typeName)
{
  xmlSetGenericErrorFunc(nullptr, xmlErrorSuppress);
  xmlDocPtr doc = xmlReadFile(std::string(filename).c_str(), nullptr, 
                               XML_PARSE_NOWARNING | XML_PARSE_NOERROR);
  if (!doc) return false;
  
  auto cleanup = [&]() { xmlFreeDoc(doc); };
  
  xmlNodePtr root = xmlDocGetRootElement(doc);
  if (!root || xmlStrcmp(root->name, BAD_CAST "VTKFile") != 0) {
    cleanup();
    return false;
  }
  
  xmlChar *type = xmlGetProp(root, BAD_CAST "type");
  bool typeMatch = type && xmlStrcmp(type, BAD_CAST typeName) == 0;
  if (type) xmlFree(type);
  if (!typeMatch) {
    cleanup();
    return false;
  }
  
  bool hasChild = false;
  for (xmlNodePtr child = root->children; child; child = child->next) {
    if (child->type == XML_ELEMENT_NODE && 
        xmlStrcmp(child->name, BAD_CAST typeName) == 0) {
      hasChild = true;
      break;
    }
  }
  
  cleanup();
  return hasChild;
}

bool validateVTUFile(std::string_view filename) {
  return validateVTKFile(filename, "UnstructuredGrid");
}

bool validateVTPFile(std::string_view filename) {
  return validateVTKFile(filename, "PolyData");
}

ErrorPredicate errorContains(std::string_view substring)
{
  return [substring](const ::precice::Error &e) -> bool {
    std::string msg(e.what());
    return msg.find(substring) != std::string::npos;
  };
}

ErrorPredicate errorMatches(std::string pattern)
{
  return [regex = std::regex(pattern)](const ::precice::Error &e) -> bool {
    return std::regex_search(e.what(), regex);
  };
}

} // namespace precice::testing
