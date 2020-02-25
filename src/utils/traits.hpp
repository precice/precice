#pragma once

#include <type_traits>

namespace precice {
namespace utils {

// Implementation of C++17 std::conjunction taken from https://en.cppreference.com/w/cpp/types/conjunction
template <class...>
struct conjunction : std::true_type {
};

template <class B1>
struct conjunction<B1> : B1 {
};

template <class B1, class... Bn>
struct conjunction<B1, Bn...>
    : std::conditional<bool(B1::value), conjunction<Bn...>, B1>::type {
};

template <class F, class... Inputs>
struct type_transform {
};

template <class P, class A1, class... A>
struct transform;

} // namespace utils
} // namespace precice
