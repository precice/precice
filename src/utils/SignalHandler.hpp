#pragma once

namespace precice {
namespace utils {

/// To be called when everything is lost and we try to terminate gracefully.
/**
 * - Writes out the event timings
 * - Deletes the address files
 */
void terminationSignalHandler(int signal);

}
}
