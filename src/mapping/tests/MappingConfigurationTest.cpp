#include <memory>
#include <string>
#include <vector>
#include "logging/Logger.hpp"
#include "mapping/Mapping.hpp"
#include "mapping/config/MappingConfiguration.hpp"
#include "mesh/SharedPointer.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"
#include "xml/XMLTag.hpp"

using namespace precice;
using namespace precice::mapping;

BOOST_AUTO_TEST_SUITE(MappingTests)
BOOST_AUTO_TEST_SUITE(Configuration)

BOOST_AUTO_TEST_CASE(Configuration)
{
  PRECICE_TEST(1_rank);

  std::string pathToTests = testing::getPathToSources() + "/mapping/tests/";
  std::string file(pathToTests + "mapping-config.xml");
  using xml::XMLTag;
  XMLTag                        tag = xml::getRootTag();
  mesh::PtrDataConfiguration    dataConfig(new mesh::DataConfiguration(tag));
  mesh::PtrMeshConfiguration    meshConfig(new mesh::MeshConfiguration(tag, dataConfig));
  mapping::MappingConfiguration mappingConfig(tag, meshConfig);
  xml::configure(tag, xml::ConfigurationContext{}, file);

  BOOST_TEST(meshConfig->meshes().size() == 3);
  BOOST_TEST(mappingConfig.mappings().size() == 3);
  BOOST_TEST(mappingConfig.mappings().at(0).fromMesh == meshConfig->meshes().at(0));
  BOOST_TEST(mappingConfig.mappings().at(0).toMesh == meshConfig->meshes().at(2));
  BOOST_TEST(mappingConfig.mappings().at(0).direction == MappingConfiguration::WRITE);

  BOOST_TEST(mappingConfig.mappings().at(1).fromMesh == meshConfig->meshes().at(2));
  BOOST_TEST(mappingConfig.mappings().at(1).toMesh == meshConfig->meshes().at(1));
  BOOST_TEST(mappingConfig.mappings().at(1).direction == MappingConfiguration::READ);

  BOOST_TEST(mappingConfig.mappings().at(2).fromMesh == meshConfig->meshes().at(1));
  BOOST_TEST(mappingConfig.mappings().at(2).toMesh == meshConfig->meshes().at(0));
  BOOST_TEST(mappingConfig.mappings().at(2).direction == MappingConfiguration::WRITE);
}

BOOST_AUTO_TEST_CASE(RBFDirectConfiguration)
{
  PRECICE_TEST(1_rank);

  std::string pathToTests = testing::getPathToSources() + "/mapping/tests/";
  std::string file(pathToTests + "mapping-rbf-direct-config.xml");
  using xml::XMLTag;
  XMLTag                        tag = xml::getRootTag();
  mesh::PtrDataConfiguration    dataConfig(new mesh::DataConfiguration(tag));
  mesh::PtrMeshConfiguration    meshConfig(new mesh::MeshConfiguration(tag, dataConfig));
  mapping::MappingConfiguration mappingConfig(tag, meshConfig);
  xml::configure(tag, xml::ConfigurationContext{}, file);

  BOOST_TEST(meshConfig->meshes().size() == 13);
  BOOST_TEST(mappingConfig.mappings().size() == 12);
  for (unsigned int i = 0; i < mappingConfig.mappings().size(); ++i) {
    BOOST_TEST(mappingConfig.mappings().at(i).mapping != nullptr);
    BOOST_TEST(mappingConfig.mappings().at(i).fromMesh == meshConfig->meshes().at(i + 1));
    BOOST_TEST(mappingConfig.mappings().at(i).toMesh == meshConfig->meshes().at(i));
    BOOST_TEST(mappingConfig.mappings().at(i).direction == MappingConfiguration::READ);
    BOOST_TEST(mappingConfig.mappings().at(i).requiresBasisFunction == true);
    BOOST_TEST(mappingConfig.mappings().at(i).configuredWithAliasTag == false);
  }
  {
    // last configured RBF
    bool solverSelection = mappingConfig.rbfConfig().solver == MappingConfiguration::RBFConfiguration::SystemSolver::GlobalDirect;
    BOOST_TEST(solverSelection);
    bool poly = mappingConfig.rbfConfig().polynomial == Polynomial::OFF;
    BOOST_TEST(poly);
    BOOST_TEST(mappingConfig.rbfConfig().deadAxis[0] == true);
    BOOST_TEST(mappingConfig.rbfConfig().deadAxis[1] == false);
    BOOST_TEST(mappingConfig.rbfConfig().deadAxis[2] == true);
    BOOST_TEST(mappingConfig.rbfConfig().solverRtol == 1e-9);
  }
}

BOOST_AUTO_TEST_CASE(RBFPUMConfiguration)
{
  PRECICE_TEST(1_rank);

  std::string pathToTests = testing::getPathToSources() + "/mapping/tests/";
  std::string file(pathToTests + "mapping-rbf-pum-direct-config.xml");
  using xml::XMLTag;
  XMLTag                        tag = xml::getRootTag();
  mesh::PtrDataConfiguration    dataConfig(new mesh::DataConfiguration(tag));
  mesh::PtrMeshConfiguration    meshConfig(new mesh::MeshConfiguration(tag, dataConfig));
  mapping::MappingConfiguration mappingConfig(tag, meshConfig);
  xml::configure(tag, xml::ConfigurationContext{}, file);

  BOOST_TEST(meshConfig->meshes().size() == 13);
  BOOST_TEST(mappingConfig.mappings().size() == 12);
  for (unsigned int i = 0; i < mappingConfig.mappings().size(); ++i) {
    BOOST_TEST(mappingConfig.mappings().at(i).mapping != nullptr);
    BOOST_TEST(mappingConfig.mappings().at(i).fromMesh == meshConfig->meshes().at(i + 1));
    BOOST_TEST(mappingConfig.mappings().at(i).toMesh == meshConfig->meshes().at(i));
    BOOST_TEST(mappingConfig.mappings().at(i).direction == MappingConfiguration::READ);
    BOOST_TEST(mappingConfig.mappings().at(i).requiresBasisFunction == true);
    BOOST_TEST(mappingConfig.mappings().at(i).configuredWithAliasTag == false);
  }
  {
    // last configured RBF
    bool solverSelection = mappingConfig.rbfConfig().solver == MappingConfiguration::RBFConfiguration::SystemSolver::PUMDirect;
    BOOST_TEST(solverSelection);
    bool poly = mappingConfig.rbfConfig().polynomial == Polynomial::OFF;
    BOOST_TEST(poly);
    BOOST_TEST(mappingConfig.rbfConfig().solverRtol == 1e-9);
    BOOST_TEST(mappingConfig.rbfConfig().verticesPerCluster == 10);
    BOOST_TEST(mappingConfig.rbfConfig().relativeOverlap == 0.4);
    BOOST_TEST(mappingConfig.rbfConfig().projectToInput == true);
  }
}

#ifndef PRECICE_NO_PETSC

BOOST_AUTO_TEST_CASE(RBFIterativeConfiguration)
{
  PRECICE_TEST(1_rank, Require::PETSc);

  std::string pathToTests = testing::getPathToSources() + "/mapping/tests/";
  std::string file(pathToTests + "mapping-rbf-iterative-config.xml");
  using xml::XMLTag;
  XMLTag                        tag = xml::getRootTag();
  mesh::PtrDataConfiguration    dataConfig(new mesh::DataConfiguration(tag));
  mesh::PtrMeshConfiguration    meshConfig(new mesh::MeshConfiguration(tag, dataConfig));
  mapping::MappingConfiguration mappingConfig(tag, meshConfig);
  xml::configure(tag, xml::ConfigurationContext{}, file);

  BOOST_TEST(meshConfig->meshes().size() == 13);
  BOOST_TEST(mappingConfig.mappings().size() == 12);
  for (unsigned int i = 0; i < mappingConfig.mappings().size(); ++i) {
    BOOST_TEST(mappingConfig.mappings().at(i).mapping != nullptr);
    BOOST_TEST(mappingConfig.mappings().at(i).fromMesh == meshConfig->meshes().at(i + 1));
    BOOST_TEST(mappingConfig.mappings().at(i).toMesh == meshConfig->meshes().at(i));
    BOOST_TEST(mappingConfig.mappings().at(i).direction == MappingConfiguration::WRITE);
    BOOST_TEST(mappingConfig.mappings().at(i).requiresBasisFunction == true);
    BOOST_TEST(mappingConfig.mappings().at(i).configuredWithAliasTag == false);
  }
  {
    // last configured RBF
    bool solverSelection = mappingConfig.rbfConfig().solver == MappingConfiguration::RBFConfiguration::SystemSolver::GlobalIterative;
    BOOST_TEST(solverSelection);
    bool poly = mappingConfig.rbfConfig().polynomial == Polynomial::OFF;
    BOOST_TEST(poly);
    BOOST_TEST(mappingConfig.rbfConfig().deadAxis[0] == true);
    BOOST_TEST(mappingConfig.rbfConfig().deadAxis[1] == false);
    BOOST_TEST(mappingConfig.rbfConfig().deadAxis[2] == true);
    BOOST_TEST(mappingConfig.rbfConfig().solverRtol == 1e-6);
  }
}
#endif

BOOST_AUTO_TEST_CASE(RBFAliasConfiguration)
{
  PRECICE_TEST(1_rank, Require::PETSc);

  std::string pathToTests = testing::getPathToSources() + "/mapping/tests/";
  std::string file(pathToTests + "mapping-rbf-alias-config.xml");
  using xml::XMLTag;
  XMLTag                        tag = xml::getRootTag();
  mesh::PtrDataConfiguration    dataConfig(new mesh::DataConfiguration(tag));
  mesh::PtrMeshConfiguration    meshConfig(new mesh::MeshConfiguration(tag, dataConfig));
  mapping::MappingConfiguration mappingConfig(tag, meshConfig);
  xml::configure(tag, xml::ConfigurationContext{}, file);

  BOOST_TEST(meshConfig->meshes().size() == 3);
  BOOST_TEST(mappingConfig.mappings().size() == 2);
  BOOST_TEST(mappingConfig.mappings().at(0).mapping != nullptr);
  BOOST_TEST(mappingConfig.mappings().at(0).fromMesh == meshConfig->meshes().at(0));
  BOOST_TEST(mappingConfig.mappings().at(0).toMesh == meshConfig->meshes().at(2));
  BOOST_TEST(mappingConfig.mappings().at(0).direction == MappingConfiguration::WRITE);
  BOOST_TEST(mappingConfig.mappings().at(0).requiresBasisFunction == true);
  BOOST_TEST(mappingConfig.mappings().at(0).configuredWithAliasTag == true);

  // The second mapping
  BOOST_TEST(mappingConfig.mappings().at(1).mapping != nullptr);
  BOOST_TEST(mappingConfig.mappings().at(1).fromMesh == meshConfig->meshes().at(2));
  BOOST_TEST(mappingConfig.mappings().at(1).toMesh == meshConfig->meshes().at(1));
  BOOST_TEST(mappingConfig.mappings().at(1).direction == MappingConfiguration::READ);
  BOOST_TEST(mappingConfig.mappings().at(1).requiresBasisFunction == true);
  BOOST_TEST(mappingConfig.mappings().at(1).configuredWithAliasTag == true);
  {
    // last configured RBF
    bool solverSelection = mappingConfig.rbfConfig().solver == MappingConfiguration::RBFConfiguration::SystemSolver::PUMDirect;
    BOOST_TEST(solverSelection);
    bool poly = mappingConfig.rbfConfig().polynomial == Polynomial::SEPARATE;
    BOOST_TEST(poly);
    BOOST_TEST(mappingConfig.rbfConfig().deadAxis[0] == false);
    BOOST_TEST(mappingConfig.rbfConfig().deadAxis[1] == false);
    BOOST_TEST(mappingConfig.rbfConfig().deadAxis[2] == false);
    BOOST_TEST(mappingConfig.rbfConfig().verticesPerCluster == 100);
    BOOST_TEST(mappingConfig.rbfConfig().relativeOverlap == 0.3);
    BOOST_TEST(mappingConfig.rbfConfig().projectToInput == true);
  }
}

#ifndef PRECICE_NO_GINKGO

#ifdef PRECICE_WITH_CUDA
// This test mostly runs the configuration path, but it does not test the actual
// setting in the Ginkgo solver class (being a cuda executor), since the mapping
// configuration does not expose this information
BOOST_AUTO_TEST_CASE(RBFGinkgoCudaConfiguration)
{
  PRECICE_TEST(1_rank);

  std::string pathToTests = testing::getPathToSources() + "/mapping/tests/";
  std::string file(pathToTests + "mapping-rbf-cuda-config.xml");
  using xml::XMLTag;
  XMLTag                        tag = xml::getRootTag();
  mesh::PtrDataConfiguration    dataConfig(new mesh::DataConfiguration(tag));
  mesh::PtrMeshConfiguration    meshConfig(new mesh::MeshConfiguration(tag, dataConfig));
  mapping::MappingConfiguration mappingConfig(tag, meshConfig);
  xml::configure(tag, xml::ConfigurationContext{}, file);

  BOOST_TEST(meshConfig->meshes().size() == 3);
  BOOST_TEST(mappingConfig.mappings().size() == 2);
  for (unsigned int i = 0; i < mappingConfig.mappings().size(); ++i) {
    BOOST_TEST(mappingConfig.mappings().at(i).mapping != nullptr);
    BOOST_TEST(mappingConfig.mappings().at(i).fromMesh == meshConfig->meshes().at(i + 1));
    BOOST_TEST(mappingConfig.mappings().at(i).toMesh == meshConfig->meshes().at(i));
    BOOST_TEST(mappingConfig.mappings().at(i).direction == MappingConfiguration::READ);
    BOOST_TEST(mappingConfig.mappings().at(i).requiresBasisFunction == true);
    BOOST_TEST(mappingConfig.mappings().at(i).configuredWithAliasTag == false);
  }
  {
    // last configured RBF
    bool solverSelection = mappingConfig.rbfConfig().solver == MappingConfiguration::RBFConfiguration::SystemSolver::GlobalIterative;
    BOOST_TEST(solverSelection);
    bool poly = mappingConfig.rbfConfig().polynomial == Polynomial::OFF;
    BOOST_TEST(poly);
    BOOST_TEST(mappingConfig.rbfConfig().deadAxis[0] == true);
    BOOST_TEST(mappingConfig.rbfConfig().deadAxis[1] == false);
    BOOST_TEST(mappingConfig.rbfConfig().deadAxis[2] == true);
    BOOST_TEST(mappingConfig.rbfConfig().solverRtol == 1e-6);
  }
}
#endif

#ifdef PRECICE_WITH_HIP
// This test mostly runs the configuration path, but it does not test the actual
// setting in the Ginkgo solver class (being a hip executor), since the mapping
// configuration does not expose this information
BOOST_AUTO_TEST_CASE(RBFGinkgoHipConfiguration)
{
  PRECICE_TEST(1_rank);

  std::string pathToTests = testing::getPathToSources() + "/mapping/tests/";
  std::string file(pathToTests + "mapping-rbf-hip-config.xml");
  using xml::XMLTag;
  XMLTag                        tag = xml::getRootTag();
  mesh::PtrDataConfiguration    dataConfig(new mesh::DataConfiguration(tag));
  mesh::PtrMeshConfiguration    meshConfig(new mesh::MeshConfiguration(tag, dataConfig));
  mapping::MappingConfiguration mappingConfig(tag, meshConfig);
  xml::configure(tag, xml::ConfigurationContext{}, file);

  BOOST_TEST(meshConfig->meshes().size() == 3);
  BOOST_TEST(mappingConfig.mappings().size() == 2);
  for (unsigned int i = 0; i < mappingConfig.mappings().size(); ++i) {
    BOOST_TEST(mappingConfig.mappings().at(i).mapping != nullptr);
    BOOST_TEST(mappingConfig.mappings().at(i).fromMesh == meshConfig->meshes().at(i + 1));
    BOOST_TEST(mappingConfig.mappings().at(i).toMesh == meshConfig->meshes().at(i));
    BOOST_TEST(mappingConfig.mappings().at(i).direction == MappingConfiguration::READ);
    BOOST_TEST(mappingConfig.mappings().at(i).requiresBasisFunction == true);
    BOOST_TEST(mappingConfig.mappings().at(i).configuredWithAliasTag == false);
  }
  {
    // last configured RBF
    bool solverSelection = mappingConfig.rbfConfig().solver == MappingConfiguration::RBFConfiguration::SystemSolver::GlobalIterative;
    BOOST_TEST(solverSelection);
    bool poly = mappingConfig.rbfConfig().polynomial == Polynomial::OFF;
    BOOST_TEST(poly);
    BOOST_TEST(mappingConfig.rbfConfig().deadAxis[0] == true);
    BOOST_TEST(mappingConfig.rbfConfig().deadAxis[1] == false);
    BOOST_TEST(mappingConfig.rbfConfig().deadAxis[2] == true);
    BOOST_TEST(mappingConfig.rbfConfig().solverRtol == 1e-6);
  }
}
#endif

#ifdef PRECICE_WITH_OPENMP
// This test mostly runs the configuration path, but it does not test the actual
// setting in the Ginkgo solver class (being an OpenMP executor), since the mapping
// configuration does not expose this information
BOOST_AUTO_TEST_CASE(RBFGinkgoOMPConfiguration)
{
  PRECICE_TEST(1_rank);

  std::string pathToTests = testing::getPathToSources() + "/mapping/tests/";
  std::string file(pathToTests + "mapping-rbf-omp-config.xml");
  using xml::XMLTag;
  XMLTag                        tag = xml::getRootTag();
  mesh::PtrDataConfiguration    dataConfig(new mesh::DataConfiguration(tag));
  mesh::PtrMeshConfiguration    meshConfig(new mesh::MeshConfiguration(tag, dataConfig));
  mapping::MappingConfiguration mappingConfig(tag, meshConfig);
  xml::configure(tag, xml::ConfigurationContext{}, file);

  BOOST_TEST(meshConfig->meshes().size() == 2);
  BOOST_TEST(mappingConfig.mappings().size() == 1);
  for (unsigned int i = 0; i < mappingConfig.mappings().size(); ++i) {
    BOOST_TEST(mappingConfig.mappings().at(i).mapping != nullptr);
    BOOST_TEST(mappingConfig.mappings().at(i).fromMesh == meshConfig->meshes().at(i + 1));
    BOOST_TEST(mappingConfig.mappings().at(i).toMesh == meshConfig->meshes().at(i));
    BOOST_TEST(mappingConfig.mappings().at(i).direction == MappingConfiguration::WRITE);
    BOOST_TEST(mappingConfig.mappings().at(i).requiresBasisFunction == true);
    BOOST_TEST(mappingConfig.mappings().at(i).configuredWithAliasTag == false);
  }
  {
    // last configured RBF
    bool solverSelection = mappingConfig.rbfConfig().solver == MappingConfiguration::RBFConfiguration::SystemSolver::GlobalIterative;
    BOOST_TEST(solverSelection);
    bool poly = mappingConfig.rbfConfig().polynomial == Polynomial::OFF;
    BOOST_TEST(poly);
    BOOST_TEST(mappingConfig.rbfConfig().deadAxis[0] == true);
    BOOST_TEST(mappingConfig.rbfConfig().deadAxis[1] == false);
    BOOST_TEST(mappingConfig.rbfConfig().deadAxis[2] == true);
    BOOST_TEST(mappingConfig.rbfConfig().solverRtol == 1e-6);
  }
}
#endif
#endif

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
