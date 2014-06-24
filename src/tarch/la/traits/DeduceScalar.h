// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_LA_TRAITS_DEDUCESCALAR_HPP_
#define _TARCH_LA_TRAITS_DEDUCESCALAR_HPP_

namespace tarch {
namespace la {

template<typename FirstScalar, typename SecondScalar>
struct DeduceScalar {};

template<> struct DeduceScalar<double,double> { typedef double Type; };
template<> struct DeduceScalar<float,float>   { typedef float Type; };
template<> struct DeduceScalar<int,int>       { typedef int Type; };

}} // namespace tarch, la

#endif /* _TARCH_LA_TRAITS_DEDUCESCALAR_HPP_ */
