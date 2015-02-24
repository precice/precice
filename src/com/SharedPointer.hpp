#include <memory>

namespace precice {
namespace com {

class Communication;
class CommunicationConfiguration;

using PtrCommunication              = std::shared_ptr<Communication>;
using PtrCommunicationConfiguration = std::shared_ptr<CommunicationConfiguration>;

}} // namespace precice, com
