#include "TXTWriterReaderTest.hpp"
#include "io/TXTWriter.hpp"
#include "io/TXTReader.hpp"
#include "utils/Parallel.hpp"
#include "math/math.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::io::tests::TXTWriterReaderTest)

namespace precice {
namespace io {
namespace tests {

logging::Logger TXTWriterReaderTest:: _log ( "precice::io::tests::TXTWriterReaderTest" );

TXTWriterReaderTest:: TXTWriterReaderTest()
:
  TestCase ( "precice::io::tests::TXTWriterReaderTest" )
{}

void TXTWriterReaderTest:: run()
{
  PRECICE_MASTER_ONLY {
    testMethod(test);
  }
}

void TXTWriterReaderTest:: test()
{
  TRACE();
  {
    Eigen::Matrix<double, 1, 2> output(1,2);
    TXTWriter txtWriter("TXTWriterReaderTest-matrix-1by2.txt");
    txtWriter.write(output);
    Eigen::Matrix<double, 1, 2> input(1,2);
    TXTReader txtReader("TXTWriterReaderTest-matrix-1by2.txt");
    txtReader.read(input);
    validateWithMessage(math::equals(output, input), input);
  }

  {
    Eigen::Matrix<double, 2, 1> output(1, 2);
    TXTWriter txtWriter("TXTWriterReaderTest-matrix-2by1.txt");
    txtWriter.write(output);
    Eigen::Matrix<double, 2, 1> input(1, 2); 
    TXTReader txtReader("TXTWriterReaderTest-matrix-2by1.txt");
    txtReader.read(input);
    validateWithMessage(math::equals(output, input), input);
  }

  {
    Eigen::Matrix<double, 3, 3> matOutput;
    matOutput << 1,2,3,4,5,6,7,8,9;
    TXTWriter txtWriter("TXTWriterReaderTest-matrix-3by3.txt");
    txtWriter.write(matOutput);
    
    Eigen::VectorXd vecOutput = Eigen::VectorXd::Constant(1, 1);
    txtWriter.write(vecOutput);
  
    Eigen::Matrix<double, 3, 3> matInput;
    TXTReader txtReader("TXTWriterReaderTest-matrix-3by3.txt");
    txtReader.read(matInput);
    validateWithMessage(math::equals(matOutput, matInput), matInput);

    Eigen::VectorXd vecInput(1);
    txtReader.read(vecInput);
    validateWithMessage(math::equals(vecOutput, vecInput), vecInput);
  }

  {
    Eigen::Vector3d output(1, 2, 3);
    TXTWriter txtWriter("TXTWriterReaderTest-vector-3.txt");
    txtWriter.write(output);

    Eigen::Vector3d input;
    TXTReader txtReader("TXTWriterReaderTest-vector-3.txt");
    txtReader.read(input);
    validateWithMessage(math::equals(output, input), input);
  }

}

}}} // namespace precice, io, tests
