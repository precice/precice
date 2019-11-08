include(CheckCXXSourceCompiles)

function(check_cpp_11_library)
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
    message(FATAL_ERROR "The standard library used by your compiler is not C++11 compliant, as it does not support time manipulation!")
  endif()
  unset(CPP11LibraryConforming)
  unset(CPP11CheckSource)
  unset(CMAKE_CXX_STANDARD)
  unset(CMAKE_CXX_STANDARD_REQUIRED)
endfunction(check_cpp_11_library)

check_cpp_11_library_support()
