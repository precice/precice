#include <Eigen/Core>
#include "io/TXTTableWriter.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(IOTests)

using namespace precice::io;

BOOST_AUTO_TEST_CASE(TXTTableWriterTest)
{
  PRECICE_TEST(1_rank);
  TXTTableWriter writer("io-TXTTableWriterTest.log");
  writer.addData("Timestep", TXTTableWriter::INT);
  writer.addData("Flowrate", TXTTableWriter::DOUBLE);
  writer.addData("Force2D", TXTTableWriter::VECTOR2D);
  writer.addData("Force3D", TXTTableWriter::VECTOR3D);

  for (int t = 0; t < 10; t++) {
    writer.writeData("Timestep", t);
    writer.writeData("Flowrate", 0.0 + (double) t);
    writer.writeData("Force2D", Eigen::Vector2d(0.0 + 2.0 * (double) t,
                                                0.0 + 2.0 * (double) t));
    writer.writeData("Force3D", Eigen::Vector3d(0.0 + 2.0 * (double) t,
                                                0.0 + 2.0 * (double) t,
                                                0.0 + 2.0 * (double) t));
  }
  writer.reset();

  writer.addData("Timestep", TXTTableWriter::INT);
  writer.addData("Importance", TXTTableWriter::INT);

  for (int t = 0; t < 10; t++) {
    writer.writeData("Timestep", t);
    writer.writeData("Importance", t * 2);
  }
  writer.close();
}

BOOST_AUTO_TEST_SUITE_END() // IOTests
