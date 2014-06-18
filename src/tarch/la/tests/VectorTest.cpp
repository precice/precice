#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#include "VectorTest.h"
#include "tarch/la/Vector.h"
#include "tarch/la/DynamicVector.h"
#include "tarch/la/WrappedVector.h"
#include "tarch/la/VectorAssign.h"
#include "tarch/la/VectorAssignList.h"
#include "tarch/la/VectorAssignRange.h"
#include "tarch/la/VectorOperations.h"
#include "tarch/la/VectorScalarOperations.h"
#include "tarch/la/VectorVectorOperations.h"
#include "tarch/la/VectorCompare.h"
#include <string>
#include <sstream>
#include <map>

#include "tarch/tests/TestCaseFactory.h"
registerTest(tarch::la::VectorTest)

namespace tarch {
namespace la {

VectorTest::VectorTest ()
:
  TestCase ("tarch::la::VectorTest")
{}

void VectorTest::run ()
{
  testMethod (testConstruction);
  testMethod (testAssignment);
  testMethod (testVectorOperations);
  testMethod (testVectorScalarOperations);
  testMethod (testVectorVectorOperations);
  testMethod (testWrappedVector);
  testMethod (testVectorCompare);

  testMethod (testVectorVectorCompare );
}


void VectorTest::testVectorVectorCompare() {
//  Vector<2,double> vector00;
//  Vector<2,double> vector01;
//  Vector<2,double> vector10;
//  Vector<2,double> vector11;
//
//  Vector<2,double> vectorM;
//
//  std::map<tarch::la::Vector<2,double> , int, tarch::la::VectorCompare<2> >  _vertexMap;
//
//  validate( _vertexMap.empty() );
//
//  assignList(vector00) = 0.0, 0.0;
//  validate( _vertexMap.find(vector00)==_vertexMap.end() );
//
//  assignList(vectorM) = 0.0, 0.3;
//  _vertexMap[vectorM] = 2;
//  validate( !_vertexMap.empty() );
//  validate( _vertexMap.size()==1 );
//  validate( _vertexMap.find(vectorM)!=_vertexMap.end() );
//  validate( _vertexMap.find(vector00)==_vertexMap.end() );
}


void VectorTest:: testConstruction ()
{
  // Test construction with uniform initialization
  Vector<1,int> vector1(1);
  validateEquals (vector1[0], 1);

  Vector<2,int> vector2(1);
  validateEquals (vector2[0], 1);
  validateEquals (vector2[1], 1);

  DynamicVector<int> dynvector(2,1);
  validateEquals (dynvector[0], 1);
  validateEquals (dynvector[1], 1);

  // Test construction with direct individual initialization
  Vector<2,int> vector3(1, 2);
  validateEquals (vector3[0], 1);
  validateEquals (vector3[1], 2);

  Vector<3,int> vector4(1, 2, 3);
  validateEquals (vector4[0], 1);
  validateEquals (vector4[1], 2);
  validateEquals (vector4[2], 3);

  // Test copy construction
  Vector<3,int> vector5(vector4);
  validateEquals (vector5[0], 1);
  validateEquals (vector5[1], 2);
  validateEquals (vector5[2], 3);

  DynamicVector<int> dynvector2(dynvector);
  validateEquals (dynvector2[0], 1);
  validateEquals (dynvector2[1], 1);

  Vector<3,int> vector6 = vector5;
  validateEquals (vector6[0], 1);
  validateEquals (vector6[1], 2);
  validateEquals (vector6[2], 3);

  // Construct vector by copying from dynamic vector
  Vector<2,int> vector7 = dynvector2;
  validateEquals (vector7[0], 1);
  validateEquals (vector7[1], 1);
}

void VectorTest::testAssignment ()
{
  // Assignment using functionality from VectorAssign.h:
  // Assign a list of scalars
  Vector<3,int> vector;
  assignList(vector) = 1, 2, 3;
  validateEquals (vector[0], 1);
  validateEquals (vector[1], 2);
  validateEquals (vector[2], 3);

  assign(vector) = 0;
  vector = 1, 2, 3;
  validateEquals (vector[0], 1);
  validateEquals (vector[1], 2);
  validateEquals (vector[2], 3);

  // Assign values of one type of vector to another type of vector
  DynamicVector<int> dynvector(3);
  dynvector = vector;
  validateEquals (dynvector[0], 1);
  validateEquals (dynvector[1], 2);
  validateEquals (dynvector[2], 3);

  assign(dynvector) = 0;
  dynvector = 1, 2, 3;
  validateEquals (dynvector[0], 1);
  validateEquals (dynvector[1], 2);
  validateEquals (dynvector[2], 3);

  // Assign a scalar to all vector
  assign(vector) = 0;
  validateEquals (vector[0], 0);
  validateEquals (vector[1], 0);
  validateEquals (vector[2], 0);
  // Assign a scalar to part of  vector
  assignRange(vector,0,2) = 0;
  validateEquals (vector[0], 0);
  validateEquals (vector[1], 0);
  validateEquals (vector[2], 0);
  // Assign a vector to part of  vector
  Vector<3,int> vector2;
  assign(vector2) = 1;
  assignRange(vector,0,2) = vector2;
  validateEquals (vector[0], 1);
  validateEquals (vector[1], 1);
  validateEquals (vector[2], 1);
  // Assignment using Vector, DynamicVector generic assignment operators:
  assign(dynvector) = 0;
  validateEquals (dynvector[0], 0);
  validateEquals (dynvector[1], 0);
  validateEquals (dynvector[2], 0);

  // Assign a scalar to part of  DynamicVector
  assignRange(dynvector,0,2) = 0;
  validateEquals (dynvector[0], 0);
  validateEquals (dynvector[1], 0);
  validateEquals (dynvector[2], 0);
  // Assign a vector to part of  DynamicVector
  DynamicVector<int> dynvector2(3);
  dynvector2 = vector2;
  assignRange(dynvector,0,2) = dynvector2;
  validateEquals (dynvector[0], 1);
  validateEquals (dynvector[1], 1);
  validateEquals (dynvector[2], 1);

  assignList(vector) = 1, 2, 3;
  assign(dynvector) = vector;
  validateEquals (dynvector[0], 1);
  validateEquals (dynvector[1], 2);
  validateEquals (dynvector[2], 3);

  assign(vector) = 1;
  validateEquals (vector[0], 1);
  validateEquals (vector[1], 1);
  validateEquals (vector[2], 1);

  vector = dynvector;
  validateEquals (vector[0], 1);
  validateEquals (vector[1], 2);
  validateEquals (vector[2], 3);

  // From fluid
  tarch::la::Vector<4,double> u(4.0);
  validateEquals(u(0),4.0);
  validateEquals(u(1),4.0);
  validateEquals(u(2),4.0);
  validateEquals(u(3),4.0);
  // not valid any more in tarch! -> compile error
//  u = 1.0, 2.0, 3.0, 4.0;
//  validateEquals(u(0),1.0);
//  validateEquals(u(1),2.0);
//  validateEquals(u(2),3.0);
//  validateEquals(u(3),4.0);

  double a = 11.0;
  double b = 12.0;
  double c = 13.0;
  double d = 14.0;
  // not valid any more in tarch! -> compile error
//  u = a, b, c, d;
//  validateEquals(u(0),11.0);
//  validateEquals(u(1),12.0);
//  validateEquals(u(2),13.0);
//  validateEquals(u(3),14.0);

  a = 21.0;
  b = 22.0;
  c = 23.0;
  d = 24.0;
  assignList(u) = a, b, c, d;
  validateEquals(u(0),21.0);
  validateEquals(u(1),22.0);
  validateEquals(u(2),23.0);
  validateEquals(u(3),24.0);


  tarch::la::DynamicVector<double> uDyn(12,5.5);
  uDyn(1*2 + 0) = 6.5;
  uDyn(2*2 + 0) = 7.5;
  uDyn(3*2 + 0) = 8.5;
  validateEquals(uDyn(0),5.5);
  validateEquals(uDyn(1),5.5);
  validateEquals(uDyn(2),6.5);
  validateEquals(uDyn(3),5.5);
  validateEquals(uDyn(4),7.5);
  validateEquals(uDyn(5),5.5);
  validateEquals(uDyn(6),8.5);
  validateEquals(uDyn(7),5.5);
  validateEquals(uDyn(8),5.5);
  validateEquals(uDyn(9),5.5);
  validateEquals(uDyn(10),5.5);
  validateEquals(uDyn(11),5.5);

  u =  uDyn(0*2 + 0),
       uDyn(1*2 + 0),
       uDyn(2*2 + 0),
       uDyn(3*2 + 0);
  validateEquals(u(0),5.5);
  validateEquals(u(1),6.5);
  validateEquals(u(2),7.5);
  validateEquals(u(3),8.5);
}

void VectorTest::testVectorOperations()
{
  Vector<3,int> vector;
  DynamicVector<double> dynvector(3);

  assignList(vector) = 1, -2, 2;
  assignList(dynvector) = 1.0, -2.0, 2.0;

  int result = sum(vector);
  validateEquals (result, 1);
  double dynresult = sum(dynvector);
  validateEquals (dynresult, 1.0);

  result = norm1(vector);
  validateEquals (result, 5.0);
  dynresult = norm1(dynvector);
  validateEquals (dynresult, 5.0);

//  result = norm2(vector);
//  validateEquals (result, 3.0);
  dynresult = norm2(dynvector);
  validateEquals (dynresult, 3.0);

  assignList(vector) = -1, 2, -2;
  assignList(dynvector) = -1.0, 2.0, -2.0;
  abs(vector, vector);
  validateEquals (vector[0], 1);
  validateEquals (vector[1], 2);
  validateEquals (vector[2], 2);
  abs(dynvector, dynvector);
  validateEquals (dynvector[0], 1.0);
  validateEquals (dynvector[1], 2.0);
  validateEquals (dynvector[2], 2.0);

  // Test indexMin, indexMax, min, and max methods
  assignList(vector) = 1, 2, 3;
  validateEquals (indexMin(vector), 0);
  validateEquals (indexMax(vector), 2);
  validateEquals (min(vector), 1);
  validateEquals (max(vector), 3);
  assignList(dynvector) = 2, -1, 2;
  validateEquals (indexMin(dynvector), 1);
  validateEquals (indexMax(dynvector), 0);
  validateEquals (min(dynvector), -1);
  validateEquals (max(dynvector), 2);
  // Test sqrt
  Vector<3,double> vector2;
  DynamicVector<double> dynvector2(3);
  assignList(vector2) = 4, 9, 16;
  assignList(dynvector2) = 4, 9, 16;
  vector2=sqrt(vector2);
  dynvector2=sqrt(dynvector2);
  validateEquals (vector2[0], 2);
  validateEquals (vector2[1], 3);
  validateEquals (vector2[2], 4);
  validateEquals (dynvector2[0], 2);
  validateEquals (dynvector2[1], 3);
  validateEquals (dynvector2[2], 4);


  Vector<3,int> vector3;
  assignList(vector3) = 1, 2, 3;
  std::ostringstream stream;
  stream << vector3;
  std::string vector3string ( stream.str() );
  validateEquals ( vector3string, std::string("1, 2, 3") );
  DynamicVector<int> dynvector3(3);
  assignList(dynvector3) = 1, 2, 3;
  std::ostringstream dynstream;
  dynstream << dynvector3;
  std::string dynvector3string ( dynstream.str() );
  validateEquals ( dynvector3string, std::string("1, 2, 3") );

  // Test sum subvectors
  {
    Vector<3,int> subvector(1, 2, 3);
    Vector<9,int> vector;
    vector = 1, 2, 3,
             4, 5, 6,
             7, 8, 9;
    sumSubvectors(vector, subvector);
    Vector<3,int> expected(12, 15, 18);
    validateEquals(subvector, expected);
  }
}

void VectorTest:: testVectorScalarOperations ()
{
  Vector<3,int> vector(0);
  DynamicVector<int> dynvector(3, 0);

  // Compound assignment operators
  vector += 1;
  dynvector += 1;
  validateEquals (vector[0], 1);
  validateEquals (vector[1], 1);
  validateEquals (vector[2], 1);
  validateEquals (dynvector[0], 1);
  validateEquals (dynvector[1], 1);
  validateEquals (dynvector[2], 1);

  vector -= 2;
  dynvector -= 2;
  validateEquals (vector[0], -1);
  validateEquals (vector[1], -1);
  validateEquals (vector[2], -1);
  validateEquals (dynvector[0], -1);
  validateEquals (dynvector[1], -1);
  validateEquals (dynvector[2], -1);

  vector *= -2;
  dynvector *= -2;
  validateEquals (vector[0], 2);
  validateEquals (vector[1], 2);
  validateEquals (vector[2], 2);
  validateEquals (dynvector[0], 2);
  validateEquals (dynvector[1], 2);
  validateEquals (dynvector[2], 2);

  vector /= 2;
  dynvector /= 2;
  validateEquals (vector[0], 1);
  validateEquals (vector[1], 1);
  validateEquals (vector[2], 1);
  validateEquals (dynvector[0], 1);
  validateEquals (dynvector[1], 1);
  validateEquals (dynvector[2], 1);

  // Arithmetic operators
  vector = vector * 2;
  dynvector = dynvector * 2;
  validateEquals (vector[0], 2);
  validateEquals (vector[1], 2);
  validateEquals (vector[2], 2);
  validateEquals (dynvector[0], 2);
  validateEquals (dynvector[1], 2);
  validateEquals (dynvector[2], 2);

  vector = 2 * vector;
  dynvector = 2 * dynvector;
  validateEquals (vector[0], 4);
  validateEquals (vector[1], 4);
  validateEquals (vector[2], 4);
  validateEquals (dynvector[0], 4);
  validateEquals (dynvector[1], 4);
  validateEquals (dynvector[2], 4);

  vector = vector / 4;
  dynvector = dynvector / 4;
  validateEquals (vector[0], 1);
  validateEquals (vector[1], 1);
  validateEquals (vector[2], 1);
  validateEquals (dynvector[0], 1);
  validateEquals (dynvector[1], 1);
  validateEquals (dynvector[2], 1);

  vector = vector + 1;
  dynvector = dynvector + 1;
  validateEquals (vector[0], 2);
  validateEquals (vector[1], 2);
  validateEquals (vector[2], 2);
  validateEquals (dynvector[0], 2);
  validateEquals (dynvector[1], 2);
  validateEquals (dynvector[2], 2);

  vector = 1 + vector;
  dynvector = 1 + dynvector;
  validateEquals (vector[0], 3);
  validateEquals (vector[1], 3);
  validateEquals (vector[2], 3);
  validateEquals (dynvector[0], 3);
  validateEquals (dynvector[1], 3);
  validateEquals (dynvector[2], 3);

  vector = vector - 1;
  dynvector = dynvector - 1;
  validateEquals (vector[0], 2);
  validateEquals (vector[1], 2);
  validateEquals (vector[2], 2);
  validateEquals (dynvector[0], 2);
  validateEquals (dynvector[1], 2);
  validateEquals (dynvector[2], 2);

  vector = 1 - vector;
  dynvector = 1 - dynvector;
  validateEquals (vector[0], -1);
  validateEquals (vector[1], -1);
  validateEquals (vector[2], -1);
  validateEquals (dynvector[0], -1);
  validateEquals (dynvector[1], -1);
  validateEquals (dynvector[2], -1);
}

void VectorTest:: testVectorVectorOperations ()
{
  Vector<3,int> vector(1);
  DynamicVector<int> dynvector(3, 1);

  // Addition assignment
  vector += vector;
  validateEquals (vector[0], 2);
  validateEquals (vector[1], 2);
  validateEquals (vector[2], 2);
  dynvector += dynvector;
  validateEquals (dynvector[0], 2);
  validateEquals (dynvector[1], 2);
  validateEquals (dynvector[2], 2);
  vector += dynvector;
  validateEquals (vector[0], 4);
  validateEquals (vector[1], 4);
  validateEquals (vector[2], 4);

  // Subtraction assignment
  vector -= vector;
  validateEquals (vector[0], 0);
  validateEquals (vector[1], 0);
  validateEquals (vector[2], 0);
  dynvector -= dynvector;
  validateEquals (dynvector[0], 0);
  validateEquals (dynvector[1], 0);
  validateEquals (dynvector[2], 0);
  assign(dynvector) = 1;
  vector -= dynvector;
  validateEquals (vector[0], -1);
  validateEquals (vector[1], -1);
  validateEquals (vector[2], -1);

  // Multiplication assignment
//  assign(vector) = 2, 4, 5;
//  assign(dynvector) = 2, 4, 5;
//  vector *= vector;
//  dynvector *= dynvector;
//  validateEquals (vector[0], 4);
//  validateEquals (vector[1], 16);
//  validateEquals (vector[2], 25);
//  validateEquals (dynvector[0], 4);
//  validateEquals (dynvector[1], 16);
//  validateEquals (dynvector[2], 25);
//  assign(vector) = 2, 4, 5;
//  assign(dynvector) = 2, 4, 5;
//  dynvector *= vector;
//  validateEquals (dynvector[0], 4);
//  validateEquals (dynvector[1], 16);
//  validateEquals (dynvector[2], 25);

  // Division assignment
//  assign(vector) = 10, 20, 30;
//  assign(dynvector) = 10, 20, 30;
//  vector /= vector;
//  dynvector /= dynvector;
//  validateEquals (vector[0], 1);
//  validateEquals (vector[1], 1);
//  validateEquals (vector[2], 1);
//  validateEquals (dynvector[0], 1);
//  validateEquals (dynvector[1], 1);
//  validateEquals (dynvector[2], 1);
//  assign(vector) = 10, 20, 30;
//  assign(dynvector) = 10, 20, 30;
//  vector /= dynvector;
//  validateEquals (vector[0], 1);
//  validateEquals (vector[1], 1);
//  validateEquals (vector[2], 1);

  // Addition
  assignList(vector) = 1, 2, 3;
  assignList(dynvector) = 3, 4, 5;
  vector = vector + vector;
  dynvector = dynvector + dynvector;
  validateEquals (vector[0], 2);
  validateEquals (vector[1], 4);
  validateEquals (vector[2], 6);
  validateEquals (dynvector[0], 6);
  validateEquals (dynvector[1], 8);
  validateEquals (dynvector[2], 10);
  vector = vector + dynvector;
  validateEquals (vector[0], 8);
  validateEquals (vector[1], 12);
  validateEquals (vector[2], 16);

  // Subtraction
  assignList(vector) = 1, 2, 3;
  assignList(dynvector) = 3, 4, 5;
  vector = vector - vector;
  dynvector = dynvector - dynvector;
  validateEquals (vector[0], 0);
  validateEquals (vector[1], 0);
  validateEquals (vector[2], 0);
  validateEquals (dynvector[0], 0);
  validateEquals (dynvector[1], 0);
  validateEquals (dynvector[2], 0);
  assignList(dynvector) = 1, 2, 3;
  vector = vector - dynvector;
  validateEquals (vector[0], -1);
  validateEquals (vector[1], -2);
  validateEquals (vector[2], -3);

  // Multiplication
//  assign(vector) = 2, 3, 4;
//  assign(dynvector) = 2, 3, 4;
//  vector = vector * vector;
//  dynvector = dynvector * dynvector;
//  validateEquals (vector[0], 4);
//  validateEquals (vector[1], 9);
//  validateEquals (vector[2], 16);
//  validateEquals (dynvector[0], 4);
//  validateEquals (dynvector[1], 9);
//  validateEquals (dynvector[2], 16);
//  assign(dynvector) = 2, 3, 4;
//  vector = vector * dynvector;
//  validateEquals (vector[0], 8);
//  validateEquals (vector[1], 27);
//  validateEquals (vector[2], 64);

  // Division
  assignList(vector) = 2, 3, 4;
  assignList(dynvector) = 2, 3, 4;
  vector = vector / vector;
  dynvector = dynvector / dynvector;
  validateEquals (vector[0], 1);
  validateEquals (vector[1], 1);
  validateEquals (vector[2], 1);
  validateEquals (dynvector[0], 1);
  validateEquals (dynvector[1], 1);
  validateEquals (dynvector[2], 1);
  assignList(vector) = 4, 8, 10;
  assignList(dynvector) = 2, 2, 2;
  vector = vector / dynvector;
  validateEquals (vector[0], 2);
  validateEquals (vector[1], 4);
  validateEquals (vector[2], 5);

  // Dot product
  assignList(vector) = 2, 3, 4;
  assignList(dynvector) = 2, 3, 4;
  int dotresult = dot(vector, vector);
  int dotdynresult = dot(dynvector, dynvector);
  validateEquals (dotresult, 29);
  validateEquals (dotdynresult, 29);
  dotresult = dot(vector, dynvector);
  validateEquals (dotresult, 29);
  dotresult = vector * vector;
  dotdynresult = dynvector * dynvector;
  validateEquals (dotresult, 29);
  validateEquals (dotdynresult, 29);
  dotresult = vector * dynvector;
  validateEquals (dotresult, 29);

  // Cross product (only works for same types)
  Vector<3,int> b(2, 3, 4);
  Vector<3,int> a(1, 2, 3);
  Vector<3,int> crossresult(0);
  cross(a, b, crossresult);
  validateEquals (crossresult[0], a[1]*b[2] - a[2]*b[1]);
  validateEquals (crossresult[1], a[2]*b[0] - a[0]*b[2]);
  validateEquals (crossresult[2], a[0]*b[1] - a[1]*b[0]);

  // Comparisons
  Vector<3,double> vec0(1.0, 2.0, 3.0);
  Vector<3,double> vec1(vec0);
  validate (equals(vec0,vec1));
  validate (not oneGreater(vec0,vec1));
  validate (oneGreaterEquals(vec0,vec1));
  validate (not allGreater(vec0,vec1));
  validate (not firstGreater(vec0,vec1));

  assignList(vec0) = 2.0, 2.0, 3.0;
  validate (not equals(vec0,vec1));
  validate (oneGreater(vec0,vec1));
  validate (oneGreaterEquals(vec0,vec1));
  validate (not allGreater(vec0,vec1));
  validate (firstGreater(vec0,vec1));

  assignList(vec0) = 2.0, 3.0, 4.0;
  validate (not equals(vec0,vec1));
  validate (oneGreater(vec0,vec1));
  validate (oneGreaterEquals(vec0,vec1));
  validate (allGreater(vec0,vec1));

  // up to here vec1=vec0
  const double tolerance = 1e-14;
  vec0(0) = vec1(0);
  vec0(1) = vec1(1);
  vec0(2) = vec1(2) + 0.99 * tolerance;
  validateWithParams3(equals(vec0,vec1,tolerance),vec0,vec1,vec0(2)-vec1(2));
  validate (not oneGreater(vec0,vec1,tolerance));
  validate (oneGreaterEquals(vec0,vec1));
  validate (not allGreater(vec0,vec1,tolerance));
  validateWithParams2(not firstGreater(vec0,vec1,tolerance), vec0, vec1);

  vec0(2) = vec1(2) + 10.0 * tolerance;
  validate (not equals(vec0,vec1,tolerance));
  validate (oneGreater(vec0,vec1,tolerance));
  validate (oneGreaterEquals(vec0,vec1));
  validate (not allGreater(vec0,vec1,tolerance));
  validateWithParams3 (firstGreater(vec0,vec1,tolerance), vec0, vec1,vec0(2)-vec1(2));

  assignList(vec0) = 1.0, 2.0, 3.0;
  vec0 += 10.0 * tolerance;
  validate (not equals(vec0,vec1,tolerance));
  validate (oneGreater(vec0,vec1,tolerance));
  validate (oneGreaterEquals(vec0,vec1));
  validate (allGreater(vec0,vec1,tolerance));
  validate (firstGreater(vec0,vec1,tolerance));

  assignList(vec0) = 1.0, 2.0, 3.0;
  vec0 -= 0.99 * tolerance;
  validate (oneGreaterEquals(vec0,vec1));

  assignList(vec0) = 1.0, 1.0, 4.0;
  validate (not firstGreater(vec0,vec1));
  assignList(vec0) = 2.0, 0.0, 0.0;
  validate (firstGreater(vec0,vec1));
  assignList(vec0) = 1.0, 2.0 + 10.0 * tolerance, 2.0;
  validate (firstGreater(vec0,vec1,tolerance));

  //Test equalsReturnIndex
  assignList(vec0) = 1.0, 2.0, 3.0 + 0.99 * tolerance;
  assignList(vec1) = 1.0, 2.0, 3.0 ;
  int i=equalsReturnIndex(vec0,vec1,tolerance);
  validateEquals (i,-1);
}

void VectorTest::testWrappedVector ()
{
  Vector<3,int> vector(0);
  int raw[3] = {1, 2, 3};
  assign(vector) = wrap<3>(raw);
  validateEquals (vector[0], 1);
  validateEquals (vector[1], 2);
  validateEquals (vector[2], 3);

  int result = vector * wrap<3>(raw);
  validateEquals (result, 14);

  wrap<3>(raw) = DynamicVector<int>(3,10);
  validateEquals (raw[0], 10);
  validateEquals (raw[1], 10);
  validateEquals (raw[2], 10);

  Vector<6,int> longVector(1);
  assignList(vector) = 2, 3, 4;
//  assign(wrap<3>(&longVector[3])) = vector;
  slice<3>(longVector,3) = vector;
  validateEquals (longVector[0], 1);
  validateEquals (longVector[1], 1);
  validateEquals (longVector[2], 1);
  validateEquals (longVector[3], 2);
  validateEquals (longVector[4], 3);
  validateEquals (longVector[5], 4);

  wrap<3,int>(raw) = slice<3>(longVector,1);
  validateEquals (raw[0], 1);
  validateEquals (raw[1], 1);
  validateEquals (raw[2], 2);

  const Vector<5,int> cvector(1);
  const Vector<5,int>& rcvector = cvector;
  assign(vector) = wrap<3,int>(&rcvector[2]);

  assignList(vector) = 3, 2, 1;
  const int * data = tarch::la::raw(vector);
  wrap<3,int>(raw) = wrap<3,int>(data);
  validateEquals (raw[0], 3);
  validateEquals (raw[1], 2);
  validateEquals (raw[2], 1);
}

void VectorTest::testVectorCompare ()
{
  typedef Vector<2,double> Vector;
  std::map<Vector,int,VectorCompare<2> > vectors;
  double eps = NUMERICAL_ZERO_DIFFERENCE;
  Vector vec0 (1.0, 2.0);
  Vector vec1 (1.0 + 10.0 * eps, 1.0);
  Vector vec2 (1.0 + 10.0 * eps, 1.5);
  Vector vec3 (2.0, 0.0);
//  vectors[vec2] = 2;
//  vectors[vec3] = 3;
//  vectors[vec0] = 0;
//  vectors[vec1] = 1;
//
//  validateEquals (vectors[vec0], 0);
//  validateEquals (vectors[vec1], 1);
//  validateEquals (vectors[vec2], 2);
//  validateEquals (vectors[vec3], 3);
}

void VectorTest::testVectorConversion()
{
  Vector<2,double> vector(1.1);
  Vector<2,int> intVector;

  intVector = integer(vector);
  validateEquals(intVector(0), 1);
  validateEquals(intVector(1), 1);

  vector = Double(intVector);
  validateEquals(vector(0), 1.0);
  validateEquals(vector(1), 1.0);
}

}} // namespace tarch, la
