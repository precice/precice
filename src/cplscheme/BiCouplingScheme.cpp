#include "BiCouplingScheme.hpp"

namespace precice {
namespace cplscheme {

BiCouplingScheme::BiCouplingScheme(
    double                        maxTime,
    int                           maxTimeWindows,
    double                        timeWindowSize,
    int                           validDigits,
    const std::string &           firstParticipant,
    const std::string &           secondParticipant,
    const std::string &           localParticipant,
    m2n::PtrM2N                   m2n,
    int                           maxIterations,
    CouplingMode                  cplMode,
    constants::TimesteppingMethod dtMethod)
    : BaseCouplingScheme(maxTime, maxTimeWindows, timeWindowSize, validDigits, firstParticipant,
                         secondParticipant, localParticipant, m2n, maxIterations, cplMode, dtMethod)
{}
} // namespace cplscheme
} // namespace precice