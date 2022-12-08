#pragma once
#ifndef PRECICE_NO_GINKGO

#include <boost/range/adaptor/indexed.hpp>
#include <boost/range/irange.hpp>
#include <ginkgo/ginkgo.hpp>
#include <ginkgo/kernels/kernel_declaration.hpp>
#include <numeric>
#include "mapping/config/MappingConfiguration.hpp"
#include "mesh/Mesh.hpp"
#include "precice/types.hpp"
#include "utils/Event.hpp"

// Declare Ginkgo Kernels
GKO_DECLARE_UNIFIED(template <typename ValueType> void create_rbf_system_matrix(
    std::shared_ptr<const DefaultExecutor> exec,
    const std::size_t n1, const std::size_t n2, const ValueType shape,
    const ValueType support_radius, ValueType *mtx, ValueType *support_points,
    ValueType *target_points, const bool add_polynomial, const unsigned int extra_dims = 0));

GKO_DECLARE_UNIFIED(template <typename ValueType> void fill_polynomial_matrix(
    std::shared_ptr<const DefaultExecutor> exec,
    const std::size_t n1, const std::size_t n2, ValueType *mtx, ValueType *x, const unsigned int dims = 4));

GKO_DECLARE_UNIFIED(template <typename ValueType> void extract_upper_triangular(
    std::shared_ptr<const DefaultExecutor> exec,
    ValueType *src, ValueType *dest,
    const std::size_t i, const std::size_t j, const std::size_t N));

GKO_REGISTER_UNIFIED_OPERATION(rbf_fill_operation, create_rbf_system_matrix);
GKO_REGISTER_UNIFIED_OPERATION(polynomial_fill_operation, fill_polynomial_matrix);
GKO_REGISTER_UNIFIED_OPERATION(tril_operation, extract_upper_triangular);

namespace precice {
namespace mapping {

template <typename RADIAL_BASIS_FUNCTION_T>
class GinkgoRadialBasisFctSolver {
public:
  GinkgoRadialBasisFctSolver() = default;

  ~GinkgoRadialBasisFctSolver() = default;

private:
};

} // namespace mapping
} // namespace precice

#endif // PRECICE_NO_GINKGO