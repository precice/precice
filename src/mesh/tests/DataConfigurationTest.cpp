// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "DataConfigurationTest.hpp"
#include "mesh/config/DataConfiguration.hpp"
#include "utils/Parallel.hpp"
#include "utils/xml/XMLTag.hpp"
#include "utils/Globals.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::mesh::tests::DataConfigurationTest)

namespace precice {
namespace mesh {
namespace tests {

logging::Logger DataConfigurationTest::
  _log("precice::mesh::tests::DataConfigurationTest");

DataConfigurationTest:: DataConfigurationTest ()
:
  TestCase ( "mesh::DataConfigurationTest" )
{}

void DataConfigurationTest:: run ()
{
  PRECICE_MASTER_ONLY {
    preciceTrace ( "run()" );
    std::string filename ( utils::Globals::getPathToSources()
                           + "/mesh/tests/data-config.xml" );
    int dim = 3;
    using utils::XMLTag;
    XMLTag tag = utils::getRootTag();
    DataConfiguration dataConfig(tag);
    dataConfig.setDimensions(dim);
    utils::configure(tag, filename);
    validateEquals(dataConfig.data().size(), 3);
    validateEquals(dataConfig.data()[0].name, std::string("vector-data"));
    validateEquals(dataConfig.data()[0].dimensions, 3);
    validateEquals(dataConfig.data()[1].name, std::string("floating-data"));
    validateEquals(dataConfig.data()[1].dimensions, 1);
    validateEquals(dataConfig.data()[2].name, std::string("second-vector-data"));
    validateEquals(dataConfig.data()[2].dimensions, 3);
  }
}

}}} // namespace precice, mesh, tests
