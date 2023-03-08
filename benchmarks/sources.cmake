#
# This file lists all benchmarks that will be compiles into precice-bench
#

target_sources(precice-bench
    PRIVATE
    benchmarks/main.cpp
    benchmarks/rbf-assembly-kernels.cpp
    benchmarks/write-block-vector-data.cpp
    )
