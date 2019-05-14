#pragma once

#include "m2n/SharedPointer.hpp"
#include <string>

namespace precice {
namespace m2n {
    
struct BoundM2N {
    void connectMasters();

    void connectSlaves();

    PtrM2N m2n;
    std::string localName;
    std::string remoteName;
    bool isRequesting;
    bool localServer;
};

} // namespace m2n
} // namespace precice
