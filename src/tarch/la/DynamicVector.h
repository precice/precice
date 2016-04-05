// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LA_DYNAMICVECTOR_H_
#define _TARCH_LA_DYNAMICVECTOR_H_

#include <vector>
#include "tarch/la/traits/IsVector.h"
#include "tarch/la/VectorAssign.h"
#include "tarch/la/VectorAssignList.h"
#include "tarch/la/VectorAssignRange.h"
#include "tarch/la/VectorOperations.h"
#include "tarch/la/VectorScalarOperations.h"
#include "tarch/la/VectorVectorOperations.h"
#include "Eigen/Core"

namespace tarch {
  namespace la {
    template<typename Scalar> class DynamicVector;
  }
}

/**
 * Run-time sized vector, based on heap-memory allocation.
 *
 * This class does not contain all methods required to operate with vectors. The
 * functionality is external and applies to all vector types. Please do not add
 * any functionality, which can be implemented externally, and hence applies to
 * all vector types.
 */
template<typename Scalar>
class tarch::la::DynamicVector
{
private:

  /**
   * Implementational basis holding values.
   */
  Scalar* _values;

  /**
   * Number of elements.
   */
  int _size;

public:

  /**
   * Default constructor for empty vector.
   */
  DynamicVector();

  /**
   * Constructs vector with size components, not initialized.
   */
  explicit DynamicVector (int size);

  /**
   * Constructs vector with size components initialized to given value.
   */
  DynamicVector (int size, const Scalar& initialValue);

  /**
   * Copy constructor.
   */
  DynamicVector (const DynamicVector<Scalar>& toCopy);

  /**
   * Copy constructor to copy from any vector type.
   *
   * The only way to accomplish this with enable-if is to specify a second
   * dummy argument with default value, which is (hopefully) optimized away.
   */
  template<typename VECTOR>
  DynamicVector (
    const VECTOR& toCopy,
    typename std::enable_if< IsVector<VECTOR>::value,void*>::type = NULL );

  /// Converts from a std::vector
  DynamicVector(const std::vector<Scalar> stdvector);

  /// Converts from a Eigen vector
  DynamicVector(const Eigen::VectorXd eigenVec);

  /**
   * Destructor, frees resources.
   */
  ~DynamicVector();

  /**
   * Assignment of vector of same type.
   */
  DynamicVector<Scalar>& operator= (const DynamicVector<Scalar>& toAssign);

  /**
   * Assignment operator for list of comma separated scalar values, that has to
   * match the number of components of the vector. Otherwise a runtime assertion
   * goes wrong.
   */
  VectorAssignList<DynamicVector<Scalar> > operator=(const Scalar& value);

//  template<typename Scalar>
//  template<typename Vector>
//    typename std::enable_if<IsVector<Vector>::value,
//    void
//  >::type DynamicVector<Scalar>::append (
//    const Vector& toAppend
//  ) {
//    assertion (toAppend.size() > 0, toAppend.size());
//    Scalar* oldValues = _values;
//    int oldSize = _size;
//    _size += toAppend.size();
//    _values = new Scalar[_size];
//    if (oldSize > 0) {
//      for (int i=0; i < oldSize; i++) {
//        _values[i] = oldValues[i];
//      }
//      delete[] oldValues;
//    }
//    for (int i=0; i < toAppend.size(); i++) {
//      _values[oldSize + i] = toAppend[i];
//    }
//  }


  /**
   * Assignment of vector of different type.
   */
  template<typename VECTOR>
    typename std::enable_if< IsVector<VECTOR>::value,
    DynamicVector<Scalar>&
  >::type operator= (const VECTOR& toAssign);

  /**
   * Returns the number of components of the vector.
   */
  int size() const;

  /**
   * Appends the given vector to this vector.
   */
  template<typename Vector>
    typename std::enable_if<IsVector<Vector>::value/* && EqualScalars<Vector,DynamicVector<Scalar> >::value*/,
    void
  >::type append (const Vector& toAppend);

  /**
   * Appends another element to this vector.
   */
  void append (const Scalar& toAppend);

  /**
   * Appends size times element toAppend.
   */
  void append (int size, const Scalar& toAppend);

  /**
   * Removes all elements from the vector.
   */
  void clear();

  /**
   * Returns const component at index.
   */
  const Scalar& operator[] (int index) const;

  /**
   * Returns component at index.
   */
  Scalar& operator[] (int index);

  /**
   * Returns const component at index.
   */
  const Scalar& operator() (int index) const;

  /**
   * Returns component at index.
   */
  Scalar& operator() (int index);

  void print() const;
  void printm(const char* filename) const;

  /**
   * @brief Converts to an std::vector
   *
   * Be advised that this conversion breaks references.
   */
  operator std::vector<Scalar>() const;

  /// Converts to an Eigen::Vector
  operator Eigen::Matrix<Scalar, Eigen::Dynamic, 1>() const;

  /// Allows to assign a Eigen::Vector to a tarch::la::DynamicVector
  DynamicVector<Scalar>& operator=(Eigen::Matrix<Scalar, Eigen::Dynamic, 1> eigenVec);

  // No more methods here? They are all generic free methods now!
};

#include "DynamicVector.cpph"

#endif /* _TARCH_LA_DYNAMICVECTOR_H_ */
