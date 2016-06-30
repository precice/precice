#ifndef _TARCH_LA_VECTORASSIGNRANGE_H_
#define _TARCH_LA_VECTORASSIGNRANGE_H_

#include "tarch/la/traits/IsVector.h"
#include "tarch/la/traits/VectorTraits.h"
#include "utils/assertion.hpp"

namespace tarch {
  namespace la {
    template<typename Vector> class VectorAssignRange;

    /**
     * Returns a wrapper object around a vector which can be assigned comma separated
     * lists.
     */
    template<typename Vector>
      typename std::enable_if< IsVector<Vector>::value,
      VectorAssignRange<Vector>
    >::type assignRange (Vector & vector,int start,int stop);
  }
}

/**
 * Wrapper around vector objects to enable comma separated list assignment.
 */
template<typename Vector>
class tarch::la::VectorAssignRange
{
private:

  Vector& _vector;
  int _start;
  int _stop;

public:

  typedef VectorTraits<Vector> Traits;

  /**
   * Constructor.
   */
  VectorAssignRange (Vector& vector,int start,int stop);

  /**
   * Destructor, asserts that the right amount of elements have been assigned.
   */
  ~VectorAssignRange();

  /**
   * Assigns the first element.
   */
  VectorAssignRange<Vector> & operator= (const typename Traits::Scalar & toAssign);
  VectorAssignRange<Vector> & operator= (const typename Traits::ThisVector & toAssign);
  /**
   * Assigns all futher elements.
   */
//  VectorAssignRange<Vector> & operator, (const typename Traits::Scalar & toAssign);
};

#include "tarch/la/VectorAssignRange.cpph"

#endif /* _TARCH_LA_VECTORASSIGNLIST_H_ */
