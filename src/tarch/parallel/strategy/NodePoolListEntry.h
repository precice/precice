// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_PARALLEL_STRATEGY_NODE_POOL_LIST_ENTRY_H_
#define _TARCH_PARALLEL_STRATEGY_NODE_POOL_LIST_ENTRY_H_

#ifdef Parallel
#include <mpi.h>
#endif
#include <string>

namespace tarch {
  namespace parallel {
    namespace strategy {
      class NodePoolListEntry;
    }
  }
}

/**
 * Node Pool List Entry Storing Status of Node
 *
 * One instance in the whole application has to keep track of the ranks
 * available, which of them are idle, and whether some of them belong to others,
 * i.e. have to communicate a lot/very fast with special other ranks. These
 * relations are represented by a node pool list entry. By this class.
 *
 * In the application's realisation, the node pool does not administer the
 * collection of ranks and their state itself. Instead, it delegates this job to
 * the node pool strategy. This strategy can be adopted to the concrete problem
 * and the concrete hardware, in particular supercomputer topology.
 *
 * @author Tobias Weinzierl
 * @version $Revision: 1.4 $
 */
class tarch::parallel::strategy::NodePoolListEntry {
  public:
    /**
     * Represents the state of the worker, i.e. whether it is idle or busy.
     */
    enum State {
      IDLE,
      WORKING
    };

  private:
    /**
     * Holds the rank of the process represented by this object.
     */
    int         _rank;

    /**
     * Holds the state of the process.
     */
    State       _state;

    /**
     * Machine name
     */
    std::string _name;

  public:
    /**
     * Construct one entry. By default this entry corresponds to an idle worker.
     */
    NodePoolListEntry( int rank, const std::string& name );

    virtual ~NodePoolListEntry();

    /**
     * Activates the node. Precondition: Node is idle. Thus, the local min level
     * is overwritten by the argument level and the state is set to working.
     */
    void activate();

    /**
     * The local rank is set to 0 and the state is switched to idle.
     */
    void deActivate();

    /**
     * @return Rank of process.
     */
    int getRank() const;

    /**
     * @return Name of the node the process is running on.
     */
    std::string getNodeName() const;

    /**
     * @return Is the node idle?
     */
    bool isIdle() const;

    /**
     * An element is smaller if and only if it is idle and the subsequent node
     * than is not idle.
     *
     * @return Object is smaller
     */
    bool operator<( const NodePoolListEntry& than ) const;

    /**
     * Two entries are equal if and only if their rank equals.
     */
    bool operator==( const NodePoolListEntry& than ) const;

    /**
     * Create string representation.
     */
    void toString(std::ostream& out) const;

    /**
     * Return string representation.
     */
    std::string toString() const;
};

std::ostream& operator<<( std::ostream& out, const tarch::parallel::strategy::NodePoolListEntry& entry );

#endif
