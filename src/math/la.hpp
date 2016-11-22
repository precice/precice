#include <Eigen/Dense>
#include "utils/assertion.hpp"

namespace precice {
namespace math {

/// Sums up the components of subvectors in vector into result.
/* Assumes vector has size k*size(result), i.e. the components of 
 * vector are a sequence of smaller vectors with length of result.
 * @param[in] vector Vector that is summed up
 * @param[in, out] result Vector which holds the sum and also defines the block sizes.
 * @todo Copy tests from tarch::la
 */
template <typename DerivedA, typename DerivedB>
void sumSubvectors (const Eigen::MatrixBase<DerivedA>& vector,
                    Eigen::MatrixBase<DerivedB>& result)
{
  int vectorSize = vector.size();
  int subvectorSize = result.size();
  assertion(vectorSize > 0);
  assertion(subvectorSize > 0);
  assertion(vectorSize % subvectorSize == 0, vectorSize, subvectorSize);

  result.setZero();
  
  // Sum up subvectors
  for (int i = 0; i < vectorSize; i+=subvectorSize) {
    result += vector.block(i, 0, subvectorSize, 1);
  }      
}

}}
