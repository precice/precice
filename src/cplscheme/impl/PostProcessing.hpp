// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_POSTPROCESSING_HPP_
#define PRECICE_CPLSCHEME_POSTPROCESSING_HPP_

#include "Eigen/Dense"
#include <map>
#include <vector>
#include "../BaseCouplingScheme.hpp"
#include "../SharedPointer.hpp"

namespace precice {
   namespace cplscheme {
      struct CouplingData;
      class BaseCouplingScheme;
   }
   namespace io {
     class TXTWriter;
     class TXTReader;
   }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace cplscheme {
namespace impl {

class PostProcessing
{
public:

	static const int NOFILTER = 0;
	static const int QR1FILTER = 1;
	static const int QR1FILTER_ABS = 2;
	static const int QR2FILTER = 3;
	static const int PODFILTER = 4;

  /**
   * @brief Map from data ID to data values.
   */
  typedef std::map<int,PtrCouplingData> DataMap;

  /**
   * @brief Destructor, empty.
   */
  virtual ~PostProcessing() {}

  virtual std::vector<int> getDataIDs() const =0;

  virtual void initialize(DataMap & cpldata) =0;

  virtual void performPostProcessing(DataMap & cpldata) =0;

  virtual void iterationsConverged(DataMap & cpldata) =0;

  /**
   * @brief sets the design specification we want to meet for the objective function,
   *        i. e., we want to solve for argmin_x ||R(x) - q||, with R(x) = H(x) - x
   *        Usually we want to solve for a fixed-point of H, thus solving for argmin_x ||R(x)||
   *        with q=0.
   */
  virtual void setDesignSpecification(Eigen::VectorXd& q) =0;

  /**
   * @brief Returns the design specification for the optimization problem.
   *        Information needed to measure the convergence.
   *        In case of manifold mapping it also returns the design specification
   *        for the surrogate model which is updated in every iteration.
   */
  virtual std::map<int, utils::DynVector> getDesignSpecification(DataMap& cplData) =0; // TODO: change to call by ref when Eigen is used.


  /**
   * @brief Sets whether the solver has to evaluate the coarse or the fine model representation
   *        steers the coupling scheme and the post processing. Only needed for multilevel based PPs.
   */
  virtual void setNextModelToEvaluate(cplscheme::BaseCouplingScheme::ModelResolution* nextModel) {};

  virtual void exportState(io::TXTWriter& writer) {}

  virtual void importState(io::TXTReader& reader) {}

  /**
   * @brief performs one optimization step of the optimization problem
   *        x_k = argmin_x||f(x_k) - q_k)
   *        with the design specification q_k and the model response f(x_k)
   */
  virtual void optimize(DataMap & cplData, Eigen::VectorXd& q)
  {
   setDesignSpecification(q);
   performPostProcessing(cplData);
  };

  virtual int getDeletedColumns() {return 0;}

};

}}} // namespace precice, cplscheme, impl

#endif /* PRECICE_CPLSCHEME_POSTPROCESSING_HPP_ */
