// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "TXTTableWriterTest.hpp"
#include "io/TXTTableWriter.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Helpers.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::io::tests::TXTTableWriterTest)

namespace precice {
namespace io {
namespace tests {

tarch::logging::Log TXTTableWriterTest:: _log ("precice::io::tests::TXTTableWriterTest");

TXTTableWriterTest:: TXTTableWriterTest ():
  TestCase ( "precice::io::tests::TXTTableWriterTest" )
{}

void TXTTableWriterTest:: run ()
{
  testMethod ( test );
}

void TXTTableWriterTest:: test ()
{
  preciceTrace ( "test()" );
  TXTTableWriter writer ( "TXTTableWriterTest-table.txt" );
  writer.addData ( "Timestep", TXTTableWriter::INT );
  writer.addData ( "Flowrate", TXTTableWriter::DOUBLE );
  writer.addData ( "Force2D", TXTTableWriter::VECTOR2D );
  writer.addData ( "Force3D", TXTTableWriter::VECTOR3D );

  for ( int t=0; t < 10; t++ ) {
    writer.writeData ( "Timestep", t );
    writer.writeData ( "Flowrate", 0.0 + (double)t );
    writer.writeData ( "Force2D", utils::Vector2D(0.0 + 2.0 * (double)t) );
    writer.writeData ( "Force3D", utils::Vector3D(0.0 + 2.0 * (double)t) );
  }
  writer.close ();
}

}}} // namespace precice, io, tests
