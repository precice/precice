#pragma once

#include <boost/system/error_code.hpp>
#include <vector>
#include "com/SharedPointer.hpp"

namespace precice {
namespace com {
class Request {

public:
  static void wait(std::vector<PtrRequest> &requests);

  virtual ~Request();

  /// Non-blocking test whether the requests has completed
  virtual bool test() = 0;

  /// Blocking wait until requests has completed
  virtual void wait() = 0;

  /** Returns the error code of the completed request.
   *
   * @precondition either wait() has been called or
   * test() has been called and returned true.
   */
  virtual boost::system::error_code errorCode() const = 0;
};
} // namespace com
} // namespace precice
