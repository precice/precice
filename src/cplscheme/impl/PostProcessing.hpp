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

  virtual void optimize(DataMap & cpldata, Eigen::VectorXd& q) {}; // change to abstract function if desigSpecification is implemented in all PPs

  virtual void iterationsConverged(DataMap & cpldata) =0;

  virtual void setDesignSpecification(Eigen::VectorXd& q) {}; // change to abstract function if desigSpecification is implemented in all PPs

  virtual void setNextModelToEvaluate(cplscheme::BaseCouplingScheme::ModelResolution& nextModel) {};

  virtual void exportState(io::TXTWriter& writer) {}

  virtual void importState(io::TXTReader& reader) {}

  virtual int getDeletedColumns() {return 0;}
};

}}} // namespace precice, cplscheme, impl

#endif /* PRECICE_CPLSCHEME_POSTPROCESSING_HPP_ */
