add_library(fmtlib-static STATIC src/format.cc src/os.cc)
target_include_directories(fmtlib-static PUBLIC "${CMAKE_CURRENT_LIST_DIR}/include")
set_target_properties(fmtlib-static
  PROPERTIES
  POSITION_INDEPENDENT_CODE ON
  CXX_STANDARD 17
  CXX_STANDARD_REQUIRED YES
  CXX_EXTENSIONS NO
  )
