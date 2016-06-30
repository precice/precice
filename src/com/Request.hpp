#ifndef PRECICE_COM_REQUEST_HPP_
#define PRECICE_COM_REQUEST_HPP_

#include <memory>
#include <vector>

namespace precice {
namespace com {
class Request {
public:
  using SharedPointer = std::shared_ptr<Request>;

public:
  static void wait(std::vector<SharedPointer>& requests);

  virtual ~Request();

  virtual bool test() = 0;

  virtual void wait() = 0;
};
}
} // namespace precice, com

#endif /* PRECICE_COM_REQUEST_HPP_ */
