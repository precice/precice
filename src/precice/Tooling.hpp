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

/** @brief Prints the configuration reference as Markdown.
 *
 * @param[inout] out the stream to write the result to
 */
void printConfigAsMD(std::ostream &out);

/** @brief Prints the configuration reference as DTD.
 *
 * @param[inout] out the stream to write the result to
 */
void printConfigAsDTD(std::ostream &out);

/** @brief Prints the configuration reference as XML with inlined help.
 *
 * @param[inout] out the stream to write the result to
 */
void printConfigAsXML(std::ostream &out);

} // namespace tooling

} // namespace precice
