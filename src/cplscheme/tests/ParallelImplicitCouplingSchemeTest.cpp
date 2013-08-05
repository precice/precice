// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ParallelImplicitCouplingSchemeTest.hpp"
#include "cplscheme/ParallelImplicitCouplingScheme.hpp"
#include "cplscheme/config/CouplingSchemeConfiguration.hpp"
#include "cplscheme/config/PostProcessingConfiguration.hpp"
#include "cplscheme/impl/ConvergenceMeasure.hpp"
#include "cplscheme/impl/AbsoluteConvergenceMeasure.hpp"
#include "cplscheme/impl/MinIterationConvergenceMeasure.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "cplscheme/Constants.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "mesh/config/MeshConfiguration.hpp"
#include "geometry/config/GeometryConfiguration.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "com/config/CommunicationConfiguration.hpp"
#include "utils/Parallel.hpp"
#include "utils/Globals.hpp"
#include "utils/xml/XMLTag.hpp"
#include "utils/Dimensions.hpp"
#include "tarch/la/Vector.h"
#include "tarch/la/WrappedVector.h"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::cplscheme::tests::ParallelImplicitCouplingSchemeTest)

namespace precice {
namespace cplscheme {
namespace tests {

using utils::Vector3D;

tarch::logging::Log ParallelImplicitCouplingSchemeTest::
  _log ( "precice::cplscheme::tests::ParallelImplicitCouplingSchemeTest" );


ParallelImplicitCouplingSchemeTest:: ParallelImplicitCouplingSchemeTest ()
:
  TestCase ( "precice::cplscheme::tests::ParallelImplicitCouplingSchemeTest" ),
  _pathToTests (),
  MY_WRITE_CHECKPOINT ( constants::actionWriteIterationCheckpoint() ),
  MY_READ_CHECKPOINT ( constants::actionReadIterationCheckpoint() )
{}

void ParallelImplicitCouplingSchemeTest:: setUp ()
{
  _pathToTests = utils::Globals::getPathToSources() + "/cplscheme/tests/";
}

void ParallelImplicitCouplingSchemeTest:: run ()
{
# ifndef PRECICE_NO_MPI
  PRECICE_MASTER_ONLY {
    testMethod(testParseConfigurationWithRelaxation);
  }
# endif // not PRECICE_NO_MPI
}

#ifndef PRECICE_NO_MPI

void ParallelImplicitCouplingSchemeTest:: testParseConfigurationWithRelaxation()
{
  preciceTrace("testParseConfigurationWithRelaxation()");
  using namespace mesh;

  std::string path(_pathToTests + "parallel-implicit-cplscheme-relax-const-config.xml");

  utils::XMLTag root = utils::getRootTag();
  PtrDataConfiguration dataConfig(new DataConfiguration(root));
  dataConfig->setDimensions(3);
  PtrMeshConfiguration meshConfig(new MeshConfiguration(root, dataConfig));
  meshConfig->setDimensions(3);
  com::PtrCommunicationConfiguration comConfig(
      new com::CommunicationConfiguration(root));
  CouplingSchemeConfiguration cplSchemeConfig(root, meshConfig, comConfig);

  utils::configure(root, path);
  validate(cplSchemeConfig._postProcConfig->getPostProcessing().get() != NULL);
  meshConfig->setMeshSubIDs();
}



#endif // not PRECICE_NO_MPI

}}} // namespace precice, cplscheme, tests
