#include <cstdint>

namespace precice::utils {

bool isMachineBigEndian()
{
  union {
    uint32_t i;
    char     c[4];
  } bint = {0x01020304};

  return bint.c[0] == 1;
}

} // namespace precice::utils
