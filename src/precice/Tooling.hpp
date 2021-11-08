#pragma once

#include <iosfwd>
#include <string>

namespace precice {

/** @brief Contains the preCICE tooling API
 *
 * The contained methods allow to query internal information of preCICE
 * without having to create a \ref precice::SolverInterface.
 *
 * @note These functions are not exposed via the bindings
 *
 * @see \ref precice::getVersionInformation which is exposed via the bindings
 */
namespace tooling {

/** The type of reference to generate
 * @see \ref precice::tooling::printConfigReference
 */
enum struct ConfigReferenceType {
  /// XML with inlined help
  XML = 0,
  /// DTD to check an XML
  DTD = 1,
  /// Markdown version used for the website
  Markdown = 2
};

/** @brief Generates a configuration reference
 *
 * @param[inout] out The stream to write the result to.
 * @param[in] reftype The type of reference to generate.
 *
 * @see \ref precice::tooling::ConfigReferenceType
 */
void printConfigReference(std::ostream &out, ConfigReferenceType reftype);

} // namespace tooling

} // namespace precice
