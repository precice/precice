include(CheckCXXSourceCompiles)

function(check_cpp_11_library_support)
  unset(CPP11LibraryConforming)
  set(CPP11CheckSource "
  #include <iostream>
  #include <iomanip>
  #include <ctime>
  int main() { std::time_t t = std::time(nullptr); std::tm tm = *std::localtime(&t); std::cout << std::put_time(&tm, \"%c %Z\"); }"
  )

  set(CMAKE_CXX_STANDARD 11)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
  check_cxx_source_compiles("${CPP11CheckSource}" CPP11LibraryConforming)
  if(NOT "${CPP11LibraryConforming}")
    message(FATAL_ERROR "
    The standard library used by your compiler ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION} is not C++11 compliant.
    It does not support time manipulation (N2071 and N2072)!
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
endfunction(check_cpp_11_library_support)

check_cpp_11_library_support()
