#include "cantorPairing.h"

namespace precice {
int
cantorPairing(int a, int b) {
  return b + (a + b) * (a + b + 1) / 2;
}
}
