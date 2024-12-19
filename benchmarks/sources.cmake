#
# This file lists all benchmarks that will be compiles into precice-bench
#

target_sources(precice-bench
    PRIVATE
    benchmarks/helper.hpp
    benchmarks/main.cpp
    benchmarks/mesh-index.cpp
    benchmarks/rbf-assembly-kernels.cpp
    benchmarks/write-data.cpp
    )
