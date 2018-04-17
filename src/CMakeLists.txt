file(GLOB sourcesAllNoMain
  "action/*.cpp"
  "action/config/*.cpp"

  "com/*.cpp"
  "com/impl/*.cpp"
  "com/config/*.cpp"

  "cplscheme/*.cpp"
  "cplscheme/impl/*.cpp"
  "cplscheme/config/*.cpp"
  
  "io/*.cpp"
  "io/impl/*.cpp"
  "io/config/*.cpp"
  
  "logging/*.cpp"
  "logging/config/*.cpp"

  "m2n/*.cpp"
  "m2n/config/*.cpp"

  "mapping/*.cpp"
  "mapping/impl/*.cpp"
  "mapping/config/*.cpp"

  "math/*.cpp"

  "mesh/*.cpp"
  "mesh/impl/*.cpp"
  "mesh/config/*.cpp"
  
  "query/*.cpp"
  
  "partition/*.cpp"
    
  "precice/*.cpp"
  "precice/impl/*.cpp"
  "precice/config/*.cpp"
  "precice/adapters/c/*.cpp"
  "precice/adapters/fortran/*.cpp"

  "utils/*.cpp"

  "xml/*.cpp"
  )

file(GLOB sourcesTests
  "*/tests/*.cpp")

set (sourcesAllNoMain ${sourcesAllNoMain} PARENT_SCOPE)
set (sourcesTests ${sourcesTests} PARENT_SCOPE) 


