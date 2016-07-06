#ifndef PRECICE_COM_COMMUNICATION_FACTORY_HPP_
#define PRECICE_COM_COMMUNICATION_FACTORY_HPP_

#include "Communication.hpp"

#include <memory>
#include <stdexcept>

namespace precice {
namespace com {
class CommunicationFactory {
public:
  using SharedPointer = std::shared_ptr<CommunicationFactory>;

public:
  /**
   * @brief Destructor.
   */
  virtual ~CommunicationFactory(){};

  virtual Communication::SharedPointer newCommunication() = 0;

  virtual std::string
  addressDirectory() {
    throw std::runtime_error("Not available!");
  }
};
}
} // namespace precice, com

#endif /* PRECICE_COM_COMMUNICATION_FACTORY_HPP_ */
