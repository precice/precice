#include "io/TXTReader.hpp"
#include "io/TXTWriter.hpp"
#include "testing/Testing.hpp"

BOOST_AUTO_TEST_SUITE(IOTests)

using namespace precice;
using namespace precice::io;

BOOST_AUTO_TEST_CASE(TXTWriterReaderTest, * testing::OnMaster())
{
  {
    Eigen::Matrix<double, 1, 2> output(1, 2);
    TXTWriter                   txtWriter("io-TXTWriterReaderTest-matrix-1by2.txt");
    txtWriter.write(output);
    Eigen::Matrix<double, 1, 2> input(1, 2);
    TXTReader                   txtReader("io-TXTWriterReaderTest-matrix-1by2.txt");
    txtReader.read(input);
    BOOST_TEST(testing::equals(output, input));
  }

  {
    Eigen::Matrix<double, 2, 1> output(1, 2);
    TXTWriter                   txtWriter("io-TXTWriterReaderTest-matrix-2by1.txt");
    txtWriter.write(output);
    Eigen::Matrix<double, 2, 1> input(1, 2);
    TXTReader                   txtReader("io-TXTWriterReaderTest-matrix-2by1.txt");
    txtReader.read(input);
    BOOST_TEST(testing::equals(output, input));
  }

  {
    Eigen::Matrix<double, 3, 3> matOutput;
    matOutput << 1, 2, 3, 4, 5, 6, 7, 8, 9;
    TXTWriter txtWriter("io-TXTWriterReaderTest-matrix-3by3.txt");
    txtWriter.write(matOutput);

    Eigen::VectorXd vecOutput = Eigen::VectorXd::Constant(1, 1);
    txtWriter.write(vecOutput);

    Eigen::Matrix<double, 3, 3> matInput;
    TXTReader                   txtReader("io-TXTWriterReaderTest-matrix-3by3.txt");
    txtReader.read(matInput);
    BOOST_TEST(testing::equals(matOutput, matInput));

    Eigen::VectorXd vecInput(1);
    txtReader.read(vecInput);
    BOOST_TEST(testing::equals(vecOutput, vecInput));
  }

  {
    Eigen::Vector3d output(1, 2, 3);
    TXTWriter       txtWriter("io-TXTWriterReaderTest-vector-3.txt");
    txtWriter.write(output);

    Eigen::Vector3d input;
    TXTReader       txtReader("io-TXTWriterReaderTest-vector-3.txt");
    txtReader.read(input);
    BOOST_TEST(testing::equals(output, input));
  }
}

BOOST_AUTO_TEST_SUITE_END() // IOTests
