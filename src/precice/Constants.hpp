#ifndef PRECICE_CONSTANTS_HPP_
#define PRECICE_CONSTANTS_HPP_

#include <string>

namespace precice {
namespace constants {

int positionInsideOfGeometry();
int positionOutsideOfGeometry();
int positionOnGeometry();

const std::string& nameConfiguration();

const std::string& dataDisplacements();
const std::string& dataForces();
const std::string& dataVelocities();

const std::string& actionWriteInitialData();
const std::string& actionWriteSimulationCheckpoint();
const std::string& actionReadSimulationCheckpoint();
const std::string& actionWriteIterationCheckpoint();
const std::string& actionReadIterationCheckpoint();
const std::string& actionPlotOutput();

int exportVTK();
int exportVRML();
int exportAll();

}} // namespace precice, constants

#endif /* PRECICE_CONSTANTS_HPP_ */
