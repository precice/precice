#include "partition/Partition.hpp"
#include "utils/MasterSlave.hpp"
#include "utils/Helpers.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Mesh.hpp"
#include "utils/MasterSlave.hpp"
#include "com/Communication.hpp"
#include <vector>
#include "utils/Globals.hpp"

namespace precice {
namespace partition {

logging::Logger Partition::_log ( "precice::partition::Partition" );

Partition:: Partition
()
:
  _mesh(),
  _fromMapping(),
  _toMapping(),
  _m2n()
{}

void Partition:: createOwnerInformation(){
  TRACE();

  if (utils::MasterSlave::_slaveMode) {
    int numberOfVertices = _mesh->vertices().size();
    utils::MasterSlave::_communication->send(numberOfVertices,0);

    if (numberOfVertices!=0) {
      std::vector<int> tags(numberOfVertices, -1);
      std::vector<int> globalIDs(numberOfVertices, -1);
      for(int i=0; i<numberOfVertices; i++){
        globalIDs[i] = _mesh->vertices()[i].getGlobalIndex();
        if(_mesh->vertices()[i].isTagged()){
          tags[i] = 1;
        }
        else{
          tags[i] = 0;
        }
      }
      DEBUG("My tags: " << tags);
      DEBUG("My global IDs: " << globalIDs);
      utils::MasterSlave::_communication->send(tags.data(),numberOfVertices,0);
      utils::MasterSlave::_communication->send(globalIDs.data(),numberOfVertices,0);
      std::vector<int> ownerVec(numberOfVertices, -1);
      utils::MasterSlave::_communication->receive(ownerVec.data(),numberOfVertices,0);
      DEBUG("My owner information: " << ownerVec);
      setOwnerInformation(ownerVec);
    }
  }


  else if (utils::MasterSlave::_masterMode) {
    std::vector<int> globalOwnerVec(_mesh->getGlobalNumberOfVertices(),0); //to temporary store which vertices already have an owner
    std::vector<std::vector<int> > slaveOwnerVecs; // the same per rank
    std::vector<std::vector<int> > slaveGlobalIDs; // global IDs per rank
    std::vector<std::vector<int> > slaveTags; // tag information per rank

    slaveOwnerVecs.resize(utils::MasterSlave::_size);
    slaveGlobalIDs.resize(utils::MasterSlave::_size);
    slaveTags.resize(utils::MasterSlave::_size);

    // fill master data

    slaveOwnerVecs[0].resize(_mesh->vertices().size());
    slaveGlobalIDs[0].resize(_mesh->vertices().size());
    slaveTags[0].resize(_mesh->vertices().size());
    for(size_t i=0; i<_mesh->vertices().size(); i++){
      slaveGlobalIDs[0][i] = _mesh->vertices()[i].getGlobalIndex();
      if(_mesh->vertices()[i].isTagged()){
        slaveTags[0][i] = 1;
      }
      else{
        slaveTags[0][i] = 0;
      }
    }

    // receive slave data
    for (int rank = 1; rank < utils::MasterSlave::_size; rank++){
      int localNumberOfVertices = -1;
      utils::MasterSlave::_communication->receive(localNumberOfVertices,rank);
      DEBUG("Rank " << rank << " has " << localNumberOfVertices << " vertices.");
      slaveOwnerVecs[rank].resize(localNumberOfVertices, 0);
      slaveTags[rank].resize(localNumberOfVertices, -1);
      slaveGlobalIDs[rank].resize(localNumberOfVertices, -1);


      if (localNumberOfVertices!=0) {
        utils::MasterSlave::_communication->receive(slaveTags[rank].data(),localNumberOfVertices,rank);
        utils::MasterSlave::_communication->receive(slaveGlobalIDs[rank].data(),localNumberOfVertices,rank);
        DEBUG("Rank " << rank << " has this tags " << slaveTags[rank]);
        DEBUG("Rank " << rank << " has this global IDs " << slaveGlobalIDs[rank]);
      }
    }

    // decide upon owners,
    int localGuess = _mesh->getGlobalNumberOfVertices() / utils::MasterSlave::_size; //guess for a decent load balancing
    //first round: every slave gets localGuess vertices
    for (int rank = 0; rank < utils::MasterSlave::_size; rank++){
      int counter = 0;
      for (size_t i=0; i < slaveOwnerVecs[rank].size(); i++) {
        if (globalOwnerVec[slaveGlobalIDs[rank][i]] == 0 && slaveTags[rank][i]==1) { // Vertex has no owner yet and rank could be owner
          slaveOwnerVecs[rank][i] = 1; // Now rank is owner
          globalOwnerVec[slaveGlobalIDs[rank][i]] = 1; //vertex now has owner
          counter++;
          if(counter==localGuess) break;
        }
      }
    }

    //second round: distribute all other vertices in a greedy way
    for (int rank = 0; rank < utils::MasterSlave::_size; rank++) {
      for(size_t i=0; i < slaveOwnerVecs[rank].size(); i++){
        if(globalOwnerVec[slaveGlobalIDs[rank][i]] == 0){
          slaveOwnerVecs[rank][i] = 1;
          globalOwnerVec[slaveGlobalIDs[rank][i]] = rank + 1;
        }
      }
    }

    // send information back to slaves
    for (int rank = 1; rank < utils::MasterSlave::_size; rank++){
      int localNumberOfVertices = slaveTags[rank].size();
      if (localNumberOfVertices!=0) {
        utils::MasterSlave::_communication->send(slaveOwnerVecs[rank].data(),localNumberOfVertices,rank);
      }
    }
    // master data
    setOwnerInformation(slaveOwnerVecs[0]);


#     ifndef NDEBUG
    for(size_t i=0;i<globalOwnerVec.size();i++){
      if(globalOwnerVec[i]==0){
        WARN( "The Vertex with global index " << i << " of mesh: " << _mesh->getName()
                       << " was completely filtered out, since it has no influence on any mapping.");
          }
    }
#     endif

  }
}

void Partition:: setOwnerInformation(const std::vector<int> &ownerVec){
  size_t i = 0;
  for ( mesh::Vertex& vertex : _mesh->vertices() ){
    assertion(i<ownerVec.size());
    assertion(ownerVec[i]!=-1);
    vertex.setOwner(ownerVec[i]==1);
    i++;
  }
}

}}
