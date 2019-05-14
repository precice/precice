#include "m2n/BoundM2N.hpp"
#include "m2n/M2N.hpp"
#include "com/Communication.hpp"

namespace precice {
namespace m2n {

void BoundM2N::connectMasters()
{
    m2n->getMasterCommunication()->prepareEstablishment();
    std::string fullLocalName = localName;
    if (localServer) fullLocalName += "Server";

    if (isRequesting){
        m2n->requestMasterConnection(remoteName, localName);
    }
    else {
        m2n->acceptMasterConnection(localName, remoteName);
    }
}

void BoundM2N::connectSlaves()
{
    if (isRequesting){
        m2n->requestSlavesConnection(remoteName, localName);
    }
    else {
        m2n->acceptSlavesConnection(localName, remoteName);
    }
    m2n->getMasterCommunication()->cleanupEstablishment();
}

} // namespace m2n
} // namespace precice
