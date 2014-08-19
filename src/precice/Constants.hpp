// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
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
