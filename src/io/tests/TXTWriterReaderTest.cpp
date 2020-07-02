#include <Eigen/Core>
#include "io/TXTReader.hpp"
#include "io/TXTWriter.hpp"
#include "logging/Logger.hpp"
#include "testing/TestContext.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(IOTests)

using namespace precice;
using namespace precice::io;

BOOST_AUTO_TEST_CASE(TXTWriterReaderTest)
{
  PRECICE_TEST(1_rank);
  {
    Eigen::Matrix<double, 1, 2> output(1, 2);
    {
      TXTWriter txtWriter("io-TXTWriterReaderTest-matrix-1by2.log");
      txtWriter.write(output);
    }
    Eigen::Matrix<double, 1, 2> input(1, 2);
    TXTReader                   txtReader("io-TXTWriterReaderTest-matrix-1by2.log");
    txtReader.read(input);
    BOOST_TEST(testing::equals(output, input));
  }

  {
    Eigen::Matrix<double, 2, 1> output(1, 2);
    {
      TXTWriter txtWriter("io-TXTWriterReaderTest-matrix-2by1.log");
      txtWriter.write(output);
    }
    Eigen::Matrix<double, 2, 1> input(1, 2);
    TXTReader                   txtReader("io-TXTWriterReaderTest-matrix-2by1.log");
    txtReader.read(input);
    BOOST_TEST(testing::equals(output, input));
  }

  {
    Eigen::Matrix<double, 3, 3> matOutput;
    matOutput << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    Eigen::VectorXd vecOutput = Eigen::VectorXd::Constant(1, 1);

    {
      TXTWriter txtWriter("io-TXTWriterReaderTest-matrix-3by3.log");
      txtWriter.write(matOutput);
      txtWriter.write(vecOutput);
    }

    Eigen::Matrix<double, 3, 3> matInput;
    TXTReader                   txtReader("io-TXTWriterReaderTest-matrix-3by3.log");
    txtReader.read(matInput);
    BOOST_TEST(testing::equals(matOutput, matInput));

    Eigen::VectorXd vecInput(1);
    txtReader.read(vecInput);
    BOOST_TEST(testing::equals(vecOutput, vecInput));
  }

  {
    Eigen::Vector3d output(1, 2, 3);
    {
      TXTWriter txtWriter("io-TXTWriterReaderTest-vector-3.log");
      txtWriter.write(output);
    }

    Eigen::Vector3d input;
    TXTReader       txtReader("io-TXTWriterReaderTest-vector-3.log");
    txtReader.read(input);
    BOOST_TEST(testing::equals(output, input));
  }
}

BOOST_AUTO_TEST_SUITE_END() // IOTests
