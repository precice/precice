#ifndef _TARCH_LA_VECTORASSIGN_H_
#define _TARCH_LA_VECTORASSIGN_H_

#include "tarch/la/traits/IsVector.h"
#include "tarch/la/traits/VectorTraits.h"
#include "utils/assertion.hpp"

namespace tarch {
  namespace la {
    template<typename Vector> class VectorAssign;

    /**
     * Returns the same vector reinterpreted to allow assignment of scalars and vectors.
     */
    template<typename Vector>
      typename std::enable_if< IsVector<Vector>::value,
      VectorAssign<Vector>&
    >::type assign (Vector& vector);
  }
}

/**
 * Reinterprets a vector allowing the assignment of scalar and vector types.
 */
template<typename Vector>
class tarch::la::VectorAssign
{
public:

   typedef VectorTraits<Vector> Traits;

   /**
    * Assigns a scalar to all components of the vector.
    */
   Vector& operator= (const typename Traits::Scalar & toAssign);

   /**
    * Assigns all components of a vector to the vector.
    */
   template<typename RVector>
     typename std::enable_if< IsVector<RVector>::value,
     Vector&
   >::type operator= (const RVector& toAssign);
};

#include "tarch/la/VectorAssign.cpph"

#endif /* _TARCH_LA_VECTORASSIGN_H_ */
