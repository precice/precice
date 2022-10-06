#pragma once

#include <boost/system/error_code.hpp>
#include <logging/Logger.hpp>

namespace precice::com {
/** Checks the error code of an async communication request.
 *
 * Will raise an error unless the error_code is success.
 *
 * @param[in] ec The error_code of the completed request
 * @param[in] _log The logger to use for the potential error message
 */
void checkErrorCode(boost::system::error_code ec, logging::Logger &_log);
} // namespace precice::com
