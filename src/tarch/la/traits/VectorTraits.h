#ifndef _TARCH_LA_TRAITS_VECTORTRAITS_H_
#define _TARCH_LA_TRAITS_VECTORTRAITS_H_

#include <vector>
#include <deque>

namespace tarch {
   namespace la {
      template<int SIZE, typename SCALAR> class Vector;
      template<typename SCALAR> class DynamicVector;
      template<int SIZE, typename SCALAR> class WrappedVector;
   }
}

namespace tarch {
namespace la {

template<typename Vector>
struct VectorTraits {};

template<int SIZE, typename SCALAR>
struct VectorTraits<Vector<SIZE,SCALAR> >
{
  typedef Vector<SIZE,SCALAR> ThisVector;
  typedef SCALAR Scalar;

  static int size (const Vector<SIZE,SCALAR>& vector)
  {
    return vector.size();
  }

  static SCALAR& elem (int index, Vector<SIZE,SCALAR>& vector )
  {
    return vector[index];
  }

  static const SCALAR& celem (int index, const Vector<SIZE,SCALAR>& vector)
  {
    return vector[index];
  }
};

template< typename SCALAR >
struct VectorTraits<DynamicVector<SCALAR> >
{
  typedef DynamicVector<SCALAR> ThisVector;
  typedef SCALAR Scalar;

  static int size ( const ThisVector& vector )
  {
    return vector.size();
  }

  static SCALAR & elem (int index, ThisVector& vector )
  {
    return vector[index];
  }

  static const SCALAR & celem ( int index, const ThisVector& vector )
  {
    return vector[index];
  }
};

template<int SIZE, typename SCALAR>
struct VectorTraits<WrappedVector<SIZE,SCALAR> >
{
   typedef WrappedVector<SIZE,SCALAR> ThisVector;
   typedef SCALAR Scalar;

   static int size ( const ThisVector& vector )
   {
      return vector.size();
   }

   static SCALAR & elem ( int index, ThisVector& vector )
   {
      return reinterpret_cast<SCALAR*>(&vector)[index];
   }

   static const SCALAR & celem (int index, const ThisVector& vector )
   {
      return reinterpret_cast<const SCALAR * const>(&vector)[index];
   }
};

template< typename SCALAR >
struct VectorTraits<std::vector<SCALAR> >
{
   typedef std::vector<SCALAR> ThisVector;
   typedef SCALAR Scalar;

   static int size (const ThisVector& vector)
   {
      return vector.size();
   }

   static SCALAR & elem (int index, ThisVector& vector)
   {
      return vector[index];
   }

   static const SCALAR & celem (int index, const ThisVector& vector )
   {
      return vector[index];
   }
};

template< typename SCALAR >
struct VectorTraits<std::deque<SCALAR> >
{
   typedef std::deque<SCALAR> ThisVector;
   typedef SCALAR Scalar;

   static int size (const ThisVector& vector)
   {
      return vector.size();
   }

   static SCALAR & elem (int index, ThisVector& vector)
   {
      return vector[index];
   }

   static const SCALAR & celem (int index, const ThisVector& vector )
   {
      return vector[index];
   }
};

}} // namespace tarch, la

#endif /* _TARCH_LA_TRAITS_VECTORTRAITS_H_ */
