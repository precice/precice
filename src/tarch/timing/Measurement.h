// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano#ifndef _PEANO_KERNEL_MULTICORE_MULTILEVELSCHEDULER_ACTION_H_
#ifndef _TARCH_TIMING_MEASUREMENT_H_
#define _TARCH_TIMING_MEASUREMENT_H_

#ifdef Parallel
#include <mpi.h>
#endif

#include <string>

#include "tarch/logging/Log.h"


namespace tarch {
  namespace timing {
    class Measurement;
  }
}


/**
 * Measurement
 *
 * The shared memory oracles have to handle real time measurements. This data
 * is always noisy, inaccurate, and has to be averaged. We originally planned
 * to move all the averaging etc. to the oracle singleton. However, that is
 * hardly possible, as the responsible oracle changes all the time, and the
 * measurements are associated to one oracle. Consequently, we introduced this
 * helper class for the individual oracle implementations.
 *
 * @author Tobias Weinzierl
 */
class tarch::timing::Measurement {
  private:
    static tarch::logging::Log _log;

    double          _accuracy;

    /**
     * Accumulated value
     *
     * The sum of all values divided by _numberOfMeasurements gives the mean
     * value.
     */
    double          _accumulatedValue;

    /**
     * To compute the standard deviation, we rely on the formula
     *
     * sigma =sqrt( E(x^2) - E(x)^2 )
     *
     * with E being the mean value.
     */
    double          _accumulatedSquares;

    /**
     * Needed to compute average value and variance.
     */
    double          _numberOfMeasurements;
    bool            _isAccurateValue;
    double          _min;
    double          _max;
    double          _minMeasurement;
    double          _maxMeasurement;

    double getMeanValue() const;
    double getMeanValueOfNextStep(double newValue) const;
  public:
    Measurement(const double& accuracy=0.0);

    /**
     * @return Averaged value (mean value) of all measurements.
     */
    double getValue() const;

    double getAccumulatedValue() const;

    double getStandardDeviation() const;

    /**
     * Is value accurate
     *
     * Whether a value is accurate depends on the last setValue() call. The
     * class internally holds the mean value of all setValue() calls. If a new
     * value is set/added, the object checks whether this additional
     * measurement modifies the mean value more than the given accuracy. If
     * this is the case, the object does not represent a valid measurement yet.
     *
     * @image Measurement.png
     *
     * It does not work if we just compare new values to the mean value. If a
     * measurement looks similar to the sketch from above (black=mean value,
     * grey area=accuracy, blue points=values), the measurement never would
     * become valid.
     *
     * I therefore use a normalised mean value
     *
     * @f$ valid = | \frac{m^{(n+1)}-m^{(n)}}{m^{(n)}} | \leq \eqs = | \frac{m^{(n+1)}}{m^{(n)}}-1  | \leq \eqs @f$
     *
     * where @f$ m^{(n)} @f$ and @f$ m^{(n+1)} @f$ are two subsequent mean
     * values.
     */
    bool isAccurateValue() const;

    /**
     * @see isAccurateValue()
     */
    void setAccuracy(const double& value);

    /**
     * Set the value. If the measurement already holds a value, this value is
     * not overwritten. Instead, the measurement accumulates all values and
     * returns the average.
     */
    void setValue(const double& value);
    int getNumberOfMeasurements() const;
    std::string toString() const;

    double max() const;
    double min() const;

    void erase();
};

#endif
