#pragma once

namespace precice::utils {

/** Ignores all inputs to prevent warnings about unused variables
 *
 * Don't worry compiers will optimize this and unneeded arguments away even at -O1.
 */
template <class... Arguments>
void ignore(Arguments &&...) {}

} // namespace precice::utils
