#include "ExportConfigurationTest.hpp"
#include "io/Export.hpp"
#include "io/config/ExportConfiguration.hpp"
#include "utils/Globals.hpp"
#include "utils/Parallel.hpp"
#include "utils/xml/XMLTag.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::io::tests::ExportConfigurationTest)

namespace precice {
namespace io {
namespace tests {


logging::Logger ExportConfigurationTest::
   _log ( "precice::io::tests::ExportConfigurationTest" );

ExportConfigurationTest:: ExportConfigurationTest ()
:
   TestCase ( "io::tests::ExportConfiguratinTest" )
{}

void ExportConfigurationTest:: run ()
{
  PRECICE_MASTER_ONLY {
    testMethod ( testConfiguration );
  }
}

void ExportConfigurationTest:: testConfiguration ()
{
  TRACE();
  using utils::XMLTag;
  XMLTag tag = utils::getRootTag();
  {
    ExportConfiguration config(tag);
    utils::configure ( tag, utils::getPathToSources() + "/io/tests/config1.xml" );
    //validate ( config.isValid() );
    validateEquals(config.exportContexts().size(), 1);
    const ExportContext& context = config.exportContexts().front();
    validateEquals ( context.type, "vtk" );
    validateEquals ( context.timestepInterval, 10 );
    //validate ( context.plotNormals );
    //validate ( context.plotNeighbors );
    validate ( context.triggerSolverPlot );
  }
  {
    tag.clear();
    ExportConfiguration config(tag);
    utils::configure ( tag, utils::getPathToSources() + "/io/tests/config2.xml" );
    //validate ( config.isValid() );
    validateEquals(config.exportContexts().size(), 1);
    const ExportContext& context = config.exportContexts().front();
    validateEquals ( context.type, "vtk" );
    validateEquals ( context.timestepInterval, 1 );
    validateEquals ( context.location, std::string("somepath") );
    //validate ( not context.plotNormals );
    //validate ( not context.plotNeighbors );
    validate ( not context.triggerSolverPlot );
  }
  {
    tag.clear();
    ExportConfiguration config(tag);
    utils::configure ( tag, utils::getPathToSources() + "/io/tests/config3.xml" );
    //validate ( config.isValid() );
    validateEquals(config.exportContexts().size(), 1);
    const ExportContext& context = config.exportContexts().front();
    validateEquals ( context.type, "vrml" );
    validateEquals ( context.timestepInterval, 1 );
    validateEquals ( context.location, std::string("") );
    //validate ( context.plotNormals );
    //validate ( not context.plotNeighbors );
    validate ( not context.triggerSolverPlot );
  }
}

}}} // namespace precice, io, tests
