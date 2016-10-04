#ifndef _TARCH_LA_VECTOR_H_
#define _TARCH_LA_VECTOR_H_

#include "tarch/la/traits/IsVector.h"
#include "tarch/la/VectorAssign.h"
#include "tarch/la/VectorAssignList.h"
#include "tarch/la/VectorAssignRange.h"
#include "tarch/la/VectorOperations.h"
#include "tarch/la/VectorScalarOperations.h"
#include "tarch/la/VectorVectorOperations.h"
#include <Eigen/Dense>

namespace tarch {
namespace la {

/**
 * Compile-time sized vector.
 *
 * This class does not contain all methods required to operate with vectors. The
 * functionality is external and applies to all vector types. Please do not add
 * any functionality, which can be implemented externally, and hence applies to
 * all vector types.
 */
template<int Size, typename Scalar>
class Vector
{
private:
  Scalar _values[Size];

public:
  Vector ();

  /**
   * Assignment operator for any vector type.
   */
  template<typename VECTOR>
    typename std::enable_if< IsVector<VECTOR>::value,
    Vector<Size,Scalar>&
  >::type operator= (const VECTOR& toAssign);

  Vector<Size,Scalar> operator= (const Eigen::Matrix<Scalar, Size, 1>& eigenVec);

  Vector<Size,Scalar> operator= (const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& eigenVec);
  
  /**
   * Assignment operator for list of comma separated scalar values, that has to
   * match the number of components of the vector. Otherwise a runtime assertion
   * goes wrong.
   */
  VectorAssignList<Vector<Size,Scalar> > operator=(const Scalar& value);

  /**
   * Copy constructor to copy from any vector type.
   *
   * The only way to accomplish this with enable-if is to specify a second
   * dummy argument with default value, which is (hopefully) optimized away.
   */
  template<typename VECTOR>
  Vector (const VECTOR& toCopy,
          typename std::enable_if< IsVector<VECTOR>::value,void*>::type = NULL);

  /**
   * Construct new vector and initialize all components with initialValue.
   */
  Vector (const Scalar& initialValue);

  /**
   * Construct new 2D vector and initialize components separately.
   *
   * Runtime assertion when dimension is not 2.
   */
  Vector (const Scalar& initialValue0,
          const Scalar& initialValue1);

  /**
   * Construct new 3D vector and initialize components separately.
   *
   * Runtime assertion when dimension is not 3.
   */
  Vector (const Scalar& initialValue0,
          const Scalar& initialValue1,
          const Scalar& initialValue2);

  /**
   * Construct new 4D vector and initialize components separately.
   *
   * Runtime assertion when dimension is not 4.
   */
  Vector (const Scalar& initialValue0,
          const Scalar& initialValue1,
          const Scalar& initialValue2,
          const Scalar& initialValue3);

  Vector (const Eigen::Matrix<Scalar, Size, 1>& eigenVec);

  Vector (const Eigen::Matrix<Scalar, Eigen::Dynamic, 1>& eigenVec);

  /**
   * Returns the number of components of the vector.
   */
  int size() const;

  /**
   * Returns read-only ref. to component of given index.
   */
  const Scalar & operator[] (int index) const;

  /**
   * Returns ref. to component of given index.
   */
  Scalar & operator[] (int index);

  /**
   * Returns read-only ref. to component of given index.
   */
  const Scalar & operator() (int index) const;

  /**
   * Returns ref. to component of given index.
   */
  Scalar & operator() (int index);

  // No more methods here? They are all generic free methods now!
};

}} // namespace tarch, la

#include "Vector.cpph"

#endif /* _TARCH_LA_VECTOR_H_ */
