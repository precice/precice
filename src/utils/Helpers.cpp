#include "utils/Globals.hpp"

#include  <boost/algorithm/string/case_conv.hpp>

namespace precice {
namespace utils {

bool isMachineBigEndian()
{
   union {
      uint32_t i;
      char c[4];
   } bint = {0x01020304};

   return bint.c[0] == 1;
}

bool convertStringToBool(std::string const & value)
{
  std::string str = value;
  boost::algorithm::to_lower(str);
  if ( str=="1" or str=="yes" or str=="true" or str=="on" )
    return true;
  
  return false;
}

}}
