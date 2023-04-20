include(CheckCXXSourceCompiles)

function(_check_cxx_11_n2071_n2072)
  unset(CXX11_N2071_N2072)
  set(CheckSourceN2071N2072 "
  #include <iostream>
  #include <iomanip>
  #include <ctime>
  int main() { std::time_t t = std::time(nullptr); std::tm tm = *std::localtime(&t); std::cout << std::put_time(&tm, \"%c %Z\"); }"
  )

  set(CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  check_cxx_source_compiles("${CheckSourceN2071N2072}" CXX11_N2071_N2072)
  if(NOT "${CXX11_N2071_N2072}")
    message(FATAL_ERROR "
    The standard library used by your compiler ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} is not C++11 compliant.
    It does not support time manipulation (N2071 and N2072).
    Please upgrade your compiler at least to:
    * GCC 5
    * LLVM 3.8
    * MSVC 19.0
    * Intel 15 in conjunction with GCC 5
    * PGI 2015 in conjunction with GCC 5
    ")
  endif()
  unset(CPP11LibraryConforming)
  unset(CPP11CheckSource)
  unset(CMAKE_CXX_STANDARD)
  unset(CMAKE_CXX_STANDARD_REQUIRED)
endfunction(_check_cxx_11_n2071_n2072)

function(_check_cxx_17_transform_reduce)
  unset(CXX17_TRANSFORM_REDUCE)
  set(CheckSourceTransformReduce "
  #include <numeric>
  #include <functional>
  #include <vector>
  int main() { std::vector<int> v{1,2,3,4,5,6}; return std::transform_reduce(v.begin(), v.end(), 0, std::plus<int>{}, [](int n){ return n * n; } ); }"
  )

  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  check_cxx_source_compiles("${CheckSourceTransformReduce}" CXX17_TRANSFORM_REDUCE)
  if(NOT "${CXX17_TRANSFORM_REDUCE}")
    message(FATAL_ERROR "
    The standard library used by your compiler ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} is not C++17 compliant.
    It does not support std::transform_reduce (P0024R2)!
    Please upgrade your compiler at least to:
    * GCC 9.3
    * LLVM 5.0
    * MSVC 19.14
    * Intel 15 in conjunction with GCC 9.2
    * PGI 2015 in conjunction with GCC 9.2
    ")
  endif()
  unset(CPP11LibraryConforming)
  unset(CPP11CheckSource)
  unset(CMAKE_CXX_STANDARD)
  unset(CMAKE_CXX_STANDARD_REQUIRED)
endfunction(_check_cxx_17_transform_reduce)

function(check_cxx_11_library_support)
  _check_cxx_11_n2071_n2072()
endfunction(check_cxx_11_library_support)

function(check_cxx_17_library_support)
  _check_cxx_17_transform_reduce()
endfunction(check_cxx_17_library_support)

check_cxx_11_library_support()
check_cxx_17_library_support()
