#include "Constants.hpp"
#include "cplscheme/Constants.hpp"
#include "io/Constants.hpp"

namespace precice {
namespace constants {

const std::string& nameConfiguration()
{
  static std::string name("precice_config.xml");
  return name;
}

const std::string& dataDisplacements()
{
  static std::string displacements("Displacements");
  return displacements;
}

const std::string& dataForces()
{
  static std::string forces("Forces");
  return forces;
}

const std::string& dataVelocities()
{
  static std::string velocities("Velocities");
  return velocities;
}

const std::string& actionWriteInitialData()
{
  return cplscheme::constants::actionWriteInitialData();
}

const std::string& actionWriteSimulationCheckpoint()
{
  static std::string writeSimulationCheckpoint("write-simulation-checkpoint");
  return writeSimulationCheckpoint;
}

const std::string& actionReadSimulationCheckpoint()
{
  static std::string readSimulationCheckpoint("read-simulation-checkpoint");
  return readSimulationCheckpoint;
}

const std::string& actionWriteIterationCheckpoint()
{
  return cplscheme::constants::actionWriteIterationCheckpoint();
}

const std::string& actionReadIterationCheckpoint()
{
  return cplscheme::constants::actionReadIterationCheckpoint();
}

const std::string& actionPlotOutput()
{
  static std::string plotOutput("plot-output");
  return plotOutput;
}

int exportVTK()
{
  return io::constants::exportVTK();
}

int exportVRML()
{
  return io::constants::exportVRML();
}

int exportAll()
{
  return io::constants::exportAll();
}

}} // namespace precice, constants
