// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_UTILS_ENABLEIF_H_
#define _TARCH_UTILS_ENABLEIF_H_

/**
 * @file
 * Enable if templates can be used to selectively activate and deactivate
 * template methods and template classes by exploiting the SFINAE (substitution-
 * failure-is-not-an-error) compiler principle (see
 * http://www.boost.org/doc/libs/1_45_0/libs/utility/enable_if.html). An example
 * on how to use enable ifs is given below.
 *
 * @code
 * template <typename Matrix>
 * struct IsColumMatrix
 * {
 *   static const bool Value = false;
 * };
 *
 * struct IsColumnMatrix<MyMatrix>
 * {
 *   static const bool Value = true;
 * };
 *
 * template <typename Matrix>
 * struct IsRowMatrix
 * {
 *   static const bool Value = false;
 * };
 *
 * template <typename Matrix>
 * EnableIf<IsColumnMatrix<Matrix> > multiply ( Matrix & A, Matrix & B );
 * {
 *   // Multiply optimized for column storage layout
 * }
 *
 * template <Matrix>
 * EnableIf< IsRowMatrix<Matrix> > multiply ( Matrix & A, Matrix & B )
 * {
 *   // Multiply optimized for row storage layout
 * }
 * @endcode
 */

namespace tarch {
namespace utils {

/**
 * Base form of enable_if, used when Enabler=true.
 */
template <bool Enabler, class ReturnType = void>
struct EnableIf {
  typedef ReturnType Type;
};

/**
 * Specialization for Enabler=false.
 */
template <class ReturnType>
struct EnableIf<false, ReturnType>
{};

/**
 * Alike enable_if_c, with nested return type only specified for Enabler=true.
 *
 * Using an enable if with invalid return type leads to an error on some compilers.
 * In order to prevent this situation, the LazyEnableIf does only extract the
 * Type information from the NestedReturnType in case the Enabler is true.
 */
template <bool Enabler, class NestedReturnType>
struct LazyEnableIf
{
  typedef typename NestedReturnType::Type Type;
};

/**
 * Specializeation for Enabler=false.
 */
template <class ReturnType>
struct LazyEnableIf<false, ReturnType>
{};

/**
 * Simple nester for usage with LazyEnableIf.
 */
template<typename TYPE>
struct LazyType
{
  typedef TYPE Type;
};

/**
 * Has a value = true, if First and Second are the same type.
 */
template<typename First, typename Second>
struct IsEqual
{
  static const int value = false;
};

/**
 * Specialization of IsEqual.
 */
template<typename Same>
struct IsEqual<Same,Same>
{
  static const int value = true;
};

}} // namespace utils, tarch

#endif /* _TARCH_UTILS_ENABLEIF_H_ */
