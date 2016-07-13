#ifndef _TARCH_MULTICORE_BOOLEAN_SEMAPHORE_OPENMP_H_
#define _TARCH_MULTICORE_BOOLEAN_SEMAPHORE_OPENMP_H_


#include <string>


namespace tarch {
  namespace multicore {
    class BooleanSemaphore;
    class Lock;
  }
}


#if defined(SharedOMP)
#include <omp.h>

class tarch::multicore::BooleanSemaphore {
  private:
    friend class tarch::multicore::Lock;

    omp_lock_t   _lock;

    /**
     * Waits until I can enter the critical section.
     */
    void enterCriticalSection();

    /**
     * Tells the semaphore that it is about to leave.
     */
    void leaveCriticalSection();

    /**
     * You may not copy a semaphore
     */
    BooleanSemaphore( const BooleanSemaphore& semaphore ) {}

    /**
     * You may not copy a semaphore
     */
    BooleanSemaphore& operator=( const BooleanSemaphore& semaphore ) {return *this;}
  public:
    BooleanSemaphore();
    ~BooleanSemaphore();

    static void sendCurrentTaskToBack(const std::string& methodTrace);
    static void continueWithTask();
};
#endif


#endif
