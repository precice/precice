#pragma once

#include <iosfwd>
#include <precice/export.h>
#include <string>

namespace precice {

/** @brief Contains the preCICE tooling API
 *
 * The contained methods allow to query internal information of preCICE
 * without having to create a \ref precice::Participant.
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
  MD = 2
};

/** @brief Generates a configuration reference
 *
 * @param[inout] out The stream to write the result to.
 * @param[in] reftype The type of reference to generate.
 *
 * @see \ref precice::tooling::ConfigReferenceType
 */
PRECICE_API void printConfigReference(std::ostream &out, ConfigReferenceType reftype);

/// @brief Checks a given configuration
PRECICE_API void checkConfiguration(const std::string &filename, const std::string &participant, int size);

} // namespace tooling

/**
 * @brief Returns information on the version of preCICE.
 *
 * Returns a semicolon-separated C-string containing:
 *
 * 1) the version of preCICE
 * 2) the revision information of preCICE
 * 3) the configuration of preCICE including MPI, PETSC, PYTHON
 */
PRECICE_API std::string getVersionInformation();

} // namespace precice
