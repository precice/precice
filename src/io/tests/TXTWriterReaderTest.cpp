// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "TXTWriterReaderTest.hpp"
#include "io/TXTWriter.hpp"
#include "io/TXTReader.hpp"
#include "tarch/la/DynamicMatrix.h"
#include "tarch/la/DynamicVector.h"
#include "utils/Globals.hpp"
#include "utils/Parallel.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::io::tests::TXTWriterReaderTest)

namespace precice {
namespace io {
namespace tests {

tarch::logging::Log TXTWriterReaderTest:: _log ( "precice::io::tests::TXTWriterReaderTest" );

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
  preciceTrace("test()");
  {
    tarch::la::DynamicMatrix<double> matrix1(1, 2);
    assignList(matrix1) = 1.0, 2.0;
    TXTWriter txtWriter("TXTWriterReaderTest-matrix-1by2.txt");
    txtWriter.write (matrix1);
  }
  {
    tarch::la::DynamicMatrix<double> matrix1(1, 2);
    tarch::la::DynamicMatrix<double> expected(1, 2);
    assignList(expected) = 1.0, 2.0;
    TXTReader txtReader("TXTWriterReaderTest-matrix-1by2.txt");
    txtReader.read(matrix1);
    validateWithMessage(tarch::la::equals(matrix1, expected), matrix1);
  }

  {
    tarch::la::DynamicMatrix<double> matrix2 (2, 1);
    assignList(matrix2) = 1.0, 2.0;
    TXTWriter txtWriter("TXTWriterReaderTest-matrix-2by1.txt");
    txtWriter.write(matrix2);
  }
  {
    tarch::la::DynamicMatrix<double> matrix2(2, 1);
    tarch::la::DynamicMatrix<double> expected(2, 1);
    assignList(expected) = 1.0, 2.0;
    TXTReader txtReader("TXTWriterReaderTest-matrix-2by1.txt");
    txtReader.read(matrix2);
    validateWithMessage(tarch::la::equals(matrix2, expected), matrix2);
  }

  {
    tarch::la::DynamicMatrix<double> matrix3 (3, 3);
    assignList(matrix3) = 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
    TXTWriter txtWriter("TXTWriterReaderTest-matrix-3by3.txt");
    txtWriter.write(matrix3);

    tarch::la::DynamicVector<double> vector1 (1);
    assignList(vector1) = 1.0;
    txtWriter.write(vector1);
  }
  {
    tarch::la::DynamicMatrix<double> matrix3 (3, 3);
    tarch::la::DynamicMatrix<double> expectedMatrix3(3, 3);
    assignList(expectedMatrix3) = 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0;
    TXTReader txtReader("TXTWriterReaderTest-matrix-3by3.txt");
    txtReader.read(matrix3);
    validateWithMessage(tarch::la::equals(matrix3, expectedMatrix3), matrix3);

    tarch::la::DynamicVector<double> vector1(1);
    tarch::la::DynamicVector<double> expectedVector1(1);
    assignList(expectedVector1) = 1.0;
    txtReader.read(vector1);
    validateWithMessage(tarch::la::equals(vector1, expectedVector1), vector1);
  }

  {
    tarch::la::DynamicVector<double> vector2(3);
    assignList(vector2) = 1.0, 2.0, 3.0;
    TXTWriter txtWriter("TXTWriterReaderTest-vector-3.txt");
    txtWriter.write(vector2);
  }
  {
    tarch::la::DynamicVector<double> vector2(3);
    tarch::la::DynamicVector<double> expected(3);
    assignList(expected) = 1.0, 2.0, 3.0;
    TXTReader txtReader("TXTWriterReaderTest-vector-3.txt");
    txtReader.read(vector2);
    validateWithMessage(tarch::la::equals(vector2, expected), vector2);
  }

}

}}} // namespace precice, io, tests
