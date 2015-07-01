#!/bin/bash

L1="$((MP_CHILD % 16))"
L2="$((L1       + 16))"

pushd 'right'
  taskset -c "${L1},${L2}" "${ATELES}" 'ateles_right.lua' &> "${configuration}-$((LEFT_NP + RIGHT_NP))-${RIGHT_NP}.log"
popd
