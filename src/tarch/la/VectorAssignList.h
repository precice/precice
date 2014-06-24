// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LA_VECTORASSIGNLIST_H_
#define _TARCH_LA_VECTORASSIGNLIST_H_

#include "tarch/la/traits/IsVector.h"
#include "tarch/la/traits/VectorTraits.h"
#include "tarch/utils/EnableIf.h"
#include "tarch/Assertions.h"

namespace tarch {
  namespace la {
    template<typename Vector> class VectorAssignList;

    /**
     * Returns a wrapper object around a vector which can be assigned comma separated
     * lists.
     */
    template<typename Vector>
      typename utils::EnableIf< IsVector<Vector>::value,
      VectorAssignList<Vector>
    >::Type assignList (Vector & vector);
  }
}

/**
 * Wrapper around vector objects to enable comma separated list assignment.
 */
template<typename Vector>
class tarch::la::VectorAssignList
{
private:

  Vector& _vector;
  int _index;

public:

  typedef VectorTraits<Vector> Traits;

  /**
   * Constructor.
   */
  VectorAssignList (Vector& vector);

  /**
   * Constructor, to start assignment not from first component (index > 0).
   */
  VectorAssignList(Vector& vector, int index);

  /**
   * Destructor, asserts that the right amount of elements have been assigned.
   */
  ~VectorAssignList();

  /**
   * Assigns the first element.
   */
  VectorAssignList<Vector> & operator= (const typename Traits::Scalar & toAssign);

  /**
   * Assigns all futher elements.
   */
  VectorAssignList<Vector> & operator, (const typename Traits::Scalar & toAssign);
};

#include "tarch/la/VectorAssignList.cpph"

#endif /* _TARCH_LA_VECTORASSIGNLIST_H_ */
