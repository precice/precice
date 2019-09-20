#pragma once

#include <Eigen/Core>
#include "utils/assertion.hpp"
#include "utils/algorithm.hpp"

namespace precice {
namespace utils {


void shiftSetFirst(Eigen::MatrixXd& A, Eigen::VectorXd& v);

void appendFront(Eigen::MatrixXd& A, Eigen::VectorXd& v);

void removeColumnFromMatrix(Eigen::MatrixXd& A, int col);

void append(Eigen::VectorXd& v, double value);

template<typename Derived1>
void append(
    Eigen::MatrixXd& A,
    const Eigen::PlainObjectBase<Derived1>& B)
{
  int n = A.rows(), m = A.cols();
  if (n <= 0 && m <= 0) {
    A = B;
  } else {
    PRECICE_ASSERT(B.rows() == n, B.rows(), A.rows());
    A.conservativeResize(n, m + B.cols());
    for(auto i = m; i < m + B.cols(); i++)
      A.col(i) = B.col(i-m);
  }
}

template<typename Derived1>
void append(
    Eigen::VectorXd& v,
    const Eigen::PlainObjectBase<Derived1>& app)
{
  int n = v.size();
  if(n <= 0){
    v = app;
  }else{
    PRECICE_ASSERT(app.cols() == 1 , app.cols());
    v.conservativeResize(n + app.size());
    for(int i = 0; i < app.size(); i++)
      v(n+i) = app(i);
  }
}

/** Maps the first at most n values of an Eigen object to a row vector
 *
 * @param[in] val the Eigen object to map from
 * @returns the mapped const row vector
 */
template<typename Derived>
auto firstN(const Eigen::PlainObjectBase<Derived>& val, unsigned n) -> const Eigen::Map<const Eigen::Matrix<typename Derived::Scalar, 1, Eigen::Dynamic>> 
{
    return {val.data(), std::min<Eigen::Index>(n,val.size())};
}


template<typename DerivedLHS, typename DerivedRHS>
bool componentWiseLess(const Eigen::PlainObjectBase<DerivedLHS>& lhs, const Eigen::PlainObjectBase<DerivedRHS>& rhs)
{
    const auto lhs_begin = lhs.data();
    const auto lhs_end   = lhs.data()+lhs.size();
    const auto rhs_begin = rhs.data();
    const auto rhs_end   = rhs.data()+rhs.size();
    
    auto mismatch = utils::mismatch(lhs_begin, lhs_end, rhs_begin, rhs_end);

    if (mismatch.first == lhs_end) {
        return true;
    }
    if (mismatch.second == rhs_end) {
        return false;
    }
    return *mismatch.first < *mismatch.second;
}

struct ComponentWiseLess {
    template<typename DerivedLHS, typename DerivedRHS>
    bool operator()(const Eigen::PlainObjectBase<DerivedLHS>& lhs, const Eigen::PlainObjectBase<DerivedRHS>& rhs) const
    {
        return componentWiseLess(lhs, rhs);
    }
};

}}
