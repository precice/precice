// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_CONVERGENCEMEASURE_HPP_
#define PRECICE_CPLSCHEME_CONVERGENCEMEASURE_HPP_

#include "cplscheme/CouplingData.hpp"
#include "utils/Dimensions.hpp"
#include "utils/Helpers.hpp"

namespace precice {
namespace cplscheme {
namespace impl {

/**
 * @brief Interface for measures checking the convergence of a series of datasets.
 *
 * A measurement involves two states of the data set: an old state and a new
 * state. Typically, the states corresponds to timestep $t$ and $t+1$. The
 * subclasses of ConvergenceMeasure define how exactly convergence is measured.
 *
 * A measure has to be used in the following way:
 * -# create the measure object (a subclass of ConvergenceMeasure)
 * -# call newMeasurementSeries() for one set of iterations
 * -# call measure() for convergence measurement
 * -# retrieve the convergence status via isConvergence()
 */
class ConvergenceMeasure
{
public:

  /**
   * @brief Destructor, empty.
   */
  virtual ~ConvergenceMeasure() {}

  /**
   * @brief To be called when a new meas. series (iteration process) starts.
   */
  virtual void newMeasurementSeries() =0;

  /**
   * @brief Performs convergence measurement.
   *
   * @param oldValues [IN] Old iterate values.
   * @param newValues [IN] New iterate values.
   */
  virtual void measure (
    const utils::DynVector& oldValues,
    const utils::DynVector& newValues ) =0;

  /**
   * @brief Returns true, if the last measurement indicates convergence.
   */
  virtual bool isConvergence() const =0;

  /**
   * @brief Adds current convergence information to output stream.
   */
  virtual std::string printState() = 0;
  
  /**
   * @brief Returns the l2-norm of the coupling residuum
   */
  virtual double getNormResidual() = 0;
};

}}} // namespace precice, cplscheme, impl

#endif /* PRECICE_CPLSCHEME_CONVERGENCEMEASURE_HPP_ */
