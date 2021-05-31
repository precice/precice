#include <Eigen/Core>
#include "utils/assertion.hpp"

namespace precice {
namespace math {

/// Sums up the components of subvectors in vector into result.
/* Assumes vector has size k*size(result), i.e. the components of 
 * vector are a sequence of smaller vectors with length of result.
 * @param[in] vector Vector that is summed up
 * @param[in, out] result Vector which holds the sum and also defines the block sizes.
 */
template <typename DerivedA, typename DerivedB>
void sumSubvectors(const Eigen::MatrixBase<DerivedA> &vector,
                   Eigen::MatrixBase<DerivedB> &      result)
{
  int vectorSize    = vector.size();
  int subvectorSize = result.size();
  PRECICE_ASSERT(vectorSize > 0);
  PRECICE_ASSERT(subvectorSize > 0);
  PRECICE_ASSERT(vectorSize % subvectorSize == 0, vectorSize, subvectorSize);

  result.setZero();

  // Sum up subvectors
  for (int i = 0; i < vectorSize; i += subvectorSize) {
    result += vector.block(i, 0, subvectorSize, 1);
  }
}

} // namespace math
} // namespace precice
