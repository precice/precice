#
# This file lists all benchmarks that will be compiles into precice-bench
#

target_sources(precice-bench
    PRIVATE
    benchmarks/bb.cpp
    benchmarks/helper.hpp
    benchmarks/main.cpp
    benchmarks/mesh-index.cpp
    benchmarks/mesh-tagging.cpp
    benchmarks/rbf-assembly-kernels.cpp
    benchmarks/write-data.cpp
    )
