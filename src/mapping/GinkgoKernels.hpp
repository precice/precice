#ifndef PRECICE_GINKGOKERNELS_HPP
#define PRECICE_GINKGOKERNELS_HPP

#include "GinkgoDefinitions.hpp"
#include "mapping/impl/BasisFunctions.hpp"


namespace precice {
namespace mapping {

namespace kernel {

template <typename EvalFunctionType>
void create_rbf_system_matrix(std::shared_ptr<const gko::Executor> exec,
                              gko::ptr_param<GinkgoMatrix> mtx, const std::array<bool, 3> activeAxis,
                              gko::ptr_param<GinkgoMatrix> supportPoints, gko::ptr_param<GinkgoMatrix> targetPoints,
                              EvalFunctionType f, ::precice::mapping::RadialBasisParameters rbf_params, bool addPolynomial, unsigned int extraDims = 0);


void fill_polynomial_matrix(std::shared_ptr<const gko::Executor> exec,
                            gko::ptr_param<GinkgoMatrix> mtx, gko::ptr_param<const GinkgoMatrix> x, const unsigned int dims);


}
}
}


#endif // PRECICE_GINKGOKERNELS_HPP
