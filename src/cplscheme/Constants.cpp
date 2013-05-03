// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "Constants.hpp"

namespace precice {
namespace cplscheme {
namespace constants {

const std::string& actionWriteIterationCheckpoint ()
{
  static std::string actionWriteIterationCheckpoint ( "write-iteration-checkpoint" );
  return actionWriteIterationCheckpoint;
}

const std::string& actionReadIterationCheckpoint ()
{
  static std::string actionReadIterationCheckpoint ( "read-iteration-checkpoint" );
  return actionReadIterationCheckpoint;
}

const std::string& actionWriteInitialData ()
{
  static std::string actionWriteInitialData ( "write-initial-data" );
  return actionWriteInitialData;
}

//const std::string WRITE_ITERATION_CHECKPOINT ( "write-iteration-checkpoint" );
//
//const std::string READ_ITERATION_CHECKPOINT ( "read-iteration-checkpoint" );
//
//const std::string WRITE_INITIAL_DATA ( "write-initial-data" );

}}} // namespace precice, cplscheme, constants
