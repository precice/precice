#!/usr/bin/env python3

import argparse
import numpy as np

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter, description='Calculates an appropriate RBF shape parameter for Gaussian basisfunctions')

# parser.add_argument('basisfunction', choices=["Gaussian"], default="Gaussian", help="Basisfunction to use")
parser.add_argument("meshwidth", type=float, help="Maximum mesh width h")
parser.add_argument("vertices", type=int, help="Vertices to include m")
parser.add_argument("--decay", type=float, default=1e-9, help="Allow basisfunction to decay to value within m vertices.")
args = parser.parse_args()

h, m, decay = args.meshwidth, args.vertices, args.decay

print("Using values:")
print("  h     =", h)
print("  m     =", m)
print("  decay =", decay)

s = np.sqrt(- np.log(decay)) / (m * h)

print("Result:")
print("  s =", s)
