// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_POSTPROCESSING_HPP_
#define PRECICE_CPLSCHEME_POSTPROCESSING_HPP_

#include <map>
#include <vector>
#include "../SharedPointer.hpp"

namespace precice {
   namespace cplscheme {
      struct CouplingData;
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

  virtual void exportState(io::TXTWriter& writer) {}

  virtual void importState(io::TXTReader& reader) {}
};

}}} // namespace precice, cplscheme, impl

#endif /* PRECICE_CPLSCHEME_POSTPROCESSING_HPP_ */
