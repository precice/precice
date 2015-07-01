#!/bin/bash

L1="$((MP_CHILD % 16))"
L2="$((L1       + 16))"

pushd 'left'
  taskset -c "${L1},${L2}" "${ATELES}" 'ateles_left.lua' &> "${configuration}-$((LEFT_NP + RIGHT_NP))-${LEFT_NP}.log"
popd
