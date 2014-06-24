#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/parallel/strategy/NodePoolListEntry.h"
#include "tarch/Assertions.h"


#include <sstream>


tarch::parallel::strategy::NodePoolListEntry::NodePoolListEntry(int rank, const std::string& name):
  _rank(rank),
  _state(WORKING),
  _name(name) {
}


tarch::parallel::strategy::NodePoolListEntry::~NodePoolListEntry() {
}


std::ostream& operator<<( std::ostream& out, const tarch::parallel::strategy::NodePoolListEntry& entry ) {
  entry.toString(out);
  return out;
}


std::string tarch::parallel::strategy::NodePoolListEntry::getNodeName() const {
  return _name;
}


std::string tarch::parallel::strategy::NodePoolListEntry::toString() const {
  std::ostringstream out;
  toString(out);
  return out.str();
}


void tarch::parallel::strategy::NodePoolListEntry::toString(std::ostream& out) const {
  out << "(rank:" << _rank;
  switch (_state) {
    case IDLE:    out << ",state:idle";     break;
    case WORKING: out << ",state:working";  break;
  }
  out << ",name:" << _name << ")";
}


bool tarch::parallel::strategy::NodePoolListEntry::operator==( const NodePoolListEntry& than ) const {
  return _rank==than._rank;
}


void tarch::parallel::strategy::NodePoolListEntry::activate() {
  assertionEquals1( _state, IDLE, toString() );
  _state = WORKING;
}


void tarch::parallel::strategy::NodePoolListEntry::deActivate() {
  _state = IDLE;
}


int tarch::parallel::strategy::NodePoolListEntry::getRank() const {
  return _rank;
}


bool tarch::parallel::strategy::NodePoolListEntry::isIdle() const {
  return _state == IDLE;
}


bool tarch::parallel::strategy::NodePoolListEntry::operator<( const tarch::parallel::strategy::NodePoolListEntry& than ) const {
  return isIdle() && !than.isIdle();
}
