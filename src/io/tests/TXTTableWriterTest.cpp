#include "TXTTableWriterTest.hpp"
#include "io/TXTTableWriter.hpp"
#include "utils/Helpers.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::io::tests::TXTTableWriterTest)

namespace precice {
namespace io {
namespace tests {

logging::Logger TXTTableWriterTest:: _log ("io::tests::TXTTableWriterTest");

TXTTableWriterTest:: TXTTableWriterTest ():
  TestCase ( "precice::io::tests::TXTTableWriterTest" )
{}

void TXTTableWriterTest:: run ()
{
  testMethod ( test );
}

void TXTTableWriterTest:: test ()
{
  TRACE();
  TXTTableWriter writer ( "TXTTableWriterTest-table.txt" );
  writer.addData ( "Timestep", TXTTableWriter::INT );
  writer.addData ( "Flowrate", TXTTableWriter::DOUBLE );
  writer.addData ( "Force2D", TXTTableWriter::VECTOR2D );
  writer.addData ( "Force3D", TXTTableWriter::VECTOR3D );

  for ( int t=0; t < 10; t++ ) {
    writer.writeData ( "Timestep", t );
    writer.writeData ( "Flowrate", 0.0 + (double)t );
    writer.writeData ( "Force2D", Eigen::Vector2d(0.0 + 2.0 * (double)t,
                                                  0.0 + 2.0 * (double)t) );
    writer.writeData ( "Force3D", Eigen::Vector3d(0.0 + 2.0 * (double)t,
                                                  0.0 + 2.0 * (double)t,
                                                  0.0 + 2.0 * (double)t) );
  }
  writer.close ();
}

}}} // namespace precice, io, tests
