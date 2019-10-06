#pragma once

#include "m2n/SharedPointer.hpp"
#include <string>

namespace precice {
namespace m2n {
    
/// An M2N between participants with a configured direction
struct BoundM2N {
    /// Prepare to establish the connection
    void prepareEstablishment();

    /// Connect the Masters of the M2N
    void connectMasters();

    /// Connect the Slaves of the M2N
    void connectSlaves();

    /// pre-connect the Slaves of the M2N
    void preConnectSlaves();

    /// Cleanup after having established the connection
    void cleanupEstablishment();

    PtrM2N m2n;
    std::string localName;
    std::string remoteName;
    bool isRequesting;
    bool localServer;
};

} // namespace m2n
} // namespace precice
