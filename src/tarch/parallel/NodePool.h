// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_PARALLEL_NODEPOOL_H_
#define _TARCH_PARALLEL_NODEPOOL_H_

#ifdef Parallel
#include <mpi.h>
#endif
#include <vector>
#include <sstream>

#include "tarch/parallel/NodePoolStrategy.h"
#include "tarch/services/Service.h"
#include "tarch/logging/Log.h"


namespace tarch {
  namespace parallel {
    class NodePool;
    class NodePoolStrategy;
  }
}

/**
 * Node Pool (Manager of Nodes)
 *
 * The node pool represents the pool of all processes available for parallel
 * computing. It holds the status of all the nodes. Hereby, the node pool of
 * the master processor is the only valid one. The other processors do not have
 * a node pool, but use the node pool as interface to access remotely the node
 * pool on rank 0. Consequently, some operations might be available on rank 0
 * only, while others are available to all ranks.
 *
 * There are basically four features implemented by this class:
 *
 * || Feature        || Operations       || For which rank || Semantics
 * |  Job management | waitForJob()      |  rank>0         |  The operation terminates as soon as the local node is assigned a job or the whole program terminates.
 * |                 | reserveFreeNode() |  all ranks      |  Get a new worker for my rank.
 * |  Hierarchy management | getMasterNodeNumber() | rank>0 | Returns the rank of the master node, i.e. of the node who has requested me as a worker.
 * |  Management     | setStrategy()     |  rank=0         |  Define which load balancing strategy to choose.
 * |                 | terminate()       |  rank=0         |  Shut down application, i.e. tell all nodes to terminate.
 *
 * The whole thing is a service on rank 0, i.e. it runs in the background of the
 * application, and regulary polls the MPI queues.
 *
 * !!! Tagging
 *
 * The node pool uses a number of tags to exchange data with the other nodes.
 * Tags are dynamically received from the node.
 *
 * || Tag || Semantics  || Message types sent on this tag
 * |  _registrationTag  | Administer overview of ranks. |  RegisterAtNodePoolMessage
 * |  _jobManagementTag | Do the rank/work assignment.  |  JobRequestMessage
 * |                    |                               |  ActivationMessage
 * |  _jobServicesTag   | Poll for new workers.         |  WorkerRequestMessage
 * |                    |                               |  NodePoolAnswerMessage
 *
 * A detailed overview of the different message types and their semantics can
 * be found no the directory's page.
 *
 * @author Tobias Weinzierl
 * @version $Revision: 1.56 $
 */
class tarch::parallel::NodePool: public tarch::services::Service {
  public:
    /**
     * Alternative return value by reserveFreeNode().
     */
    static const int NoFreeNodesMessage;

    typedef int JobRequestMessageAnswer;

    struct JobRequestMessageAnswerValues {
      /**
       * If a job request message, i.e. waitForJob() is answered with a value
       * greater or equal to zero, this means that there is a new master for
       * this job. The next thing the node should do then is to receive the
       * corresponding fork message, i.e. to inform about which data to handle.
       */
      static const int NewMaster;

      /**
       * If a job sends a job request message, i.e. calls waitForJob(), and
       * Terminate is sent back, the job should shut down.
       */
      static const int Terminate;

      /**
       * If a job sends a job request message, i.e. calls waitForJob(), and
       * HandleLocalProblem is sent back, the job should concentrate on some
       * local stuff. As soon as it has finished its local stuff, it should ask
       * for a new job.
       */
      static const int HandleLocalProblem;

      /**
       * If a job sends a job request message, i.e. calls waitForJob(), and
       * ExecuteGlobalTask is sent back, the job can be used to execute a task
       * on all processes. For example a new MPI communicator could be created
       * (which requires that all processes of an existing communicator
       * participate)
       */
      static const int ExecuteGlobalTask;
    };


  private:
    /**
     * Logging device.
     */
    static tarch::logging::Log _log;

    /**
     * Holds the number of the master process that has triggered this node to
     * run.
     */
    int _masterNode;

    int _registrationTag;
    int _jobManagementTag;
    int _jobServicesTag;

    #ifdef Asserts
    bool _isInitialised;
    #endif

    /**
     * In the beginning every instance is alive. If the instance is receiving
     * a termination message (client) or terminate() (server) is called, the
     * instance becomes not alive. Afterwards one is not allowed to ask for
     * further jobs or additional subclients.
     */
    bool _isAlive;

    bool _answerJobRequestsWithHandleLocalProblem;

    /**
     * Send ExecuteGlobalTask for all job requests
     */
    bool _answerJobRequestsWithExecuteGlobalTask;

    /**
     * Holds the strategy to use to answer the worker queries.
     */
    NodePoolStrategy* _strategy;

    /**
     * The operation reserves the number of a free process, if there's one.
     * Afterwards this free process is sent an activation message and the number
     * of the new node is returned.
     *
     * At the first glance this operation is easier compared to
     * reserveFreeNodeForClient(). Nevertheless, one has to consider a
     * sophisticated side effect: The client to be actived might have sent a
     * job request already. Thus, if the replyToMessages() operation receives
     * a job request from a node set active, it has to cancel this message
     * without handling it.
     */
    int reserveFreeNodeForServer();

    /**
     * Reserve Worker for any Process not Being Node Pool Server
     *
     * To reserve a free node for the client the operation sents a ressource
     * allocation message to the server. Afterwards it calls a blocking receive
     * operation to receive the number of a new client or a NoFreeNodesMessage.
     * This number if returned if it is a valid number, otherwise the status
     * of the node is set #not alive# and the flag NoFreeNodesMessage is
     * returned.
     *
     * If fair scheduling is switched on, it usually takes a certain time until
     * the node pool server sends back the answer. For this time the client
     * could sleep and allow other tasks to work. This is especially important
     * if a node is overloaded. Yet, if no worker is available, the node pool
     * server answers immediately and there's no reason to sleep for the asking
     * worker. Thus, I introduced the _hasReceivedNoWorkerAvailable flag passed
     * the operation. The following table gives the transitions.
     *
|| _hasReceivedNoWorkerAvailable (in) || received message    || sleep || _hasReceivedNoWorkerAvailable (out)
|  false                              |  new (sub-)worker    |  yes   |  false
|  false                              |  no worker available |  yes   |  true
|  true                               |  new (sub-)worker    |  no    |  false
|  true                               |  no worker available |  no    |  true
     *
     * To make the startup fast, the initial value of this flag should be true.
     */
    int reserveFreeNodeForClient();

    /**
     * Empty the receive buffer of the process. Is called by the constructor
     * to ensure system integrity.
     */
    void emptyReceiveBuffers();

    /**
     * Operation reserves the first free node and returns either the rank of
     * this node or NoFreeNodesMessage if there isn't one. Is invoked only on
     * the master process.
     *
     * If it was a valid number, the corresponding entry has been modified
     * before the operation terminated. Used only by server.
     */
    int getFreeNode(int forMaster);

    /**
     * One is not allowed to create a copy of a node pool.
     */
    NodePool( const NodePool& pool );

    /**
     * One is not allowed to create a copy of a node pool.
     */
    NodePool& operator=( const NodePool& pool );

    /**
     * Constructor for the root node who wants to become the job server.
     */
    NodePool();

    void replyToRegistrationMessages();
    void replyToJobRequestMessages();
    void replyToWorkerRequestMessages();

    void emptyRegisterMessageReceiveBuffer();
    void emptyJobRequestMessageBuffer();
    void emptyWorkerRequestMessageBuffer();

    /**
     * Only the master is allowed to call this operation.
     */
    void waitUntilAllTerminationMessagesHaveBeenSent();
  public:
    /**
     * Set service up.
     *
     * Should be called once and only once. Its counterpart is shutdown().
     */
    void init();

    /**
     * The node pool is a singleton. Thus, one has to use getInstance() to get
     * an instance.
     *
     * At the very beginning of the simulation, you should set a strategy on
     * the master before you do anything else in the code. See operation
     * setStrategy().
     */
    static NodePool& getInstance();

    /**
     * Destructor. Note that there is a warning if the destructor is called on
     * an alive session. Ensure it has received a termination message before to
     * avoid this warning.
     */
    virtual ~NodePool();

    /**
     * Wait for a New Job
     *
     * This operation blocks the program execution until a new job has arrived.
     * It might happen, that this operation receives a termination message
     * instead of a new job.
     *
     * As waiting for a new job might last pretty long, the software should not
     * trigger a deadlock time out. Thus, the code does not use the code
     * generated by DaStGen but invokes MPI_Send directly.
     *
     * @return What do do next
     */
    JobRequestMessageAnswer waitForJob();

    /**
     * This is a server specific operation: All processes waiting for a job are
     * told to terminate. So an internal flag is set and all new queries are
     * answered this way. Furthermore all processes marked as not working are
     * sent a termination message immediately.
     */
    void terminate();

    /**
     * Asks for a free node. The method returns the number of the node reserved
     * or -1 if there's no node available. level has to hold the level of the
     * element one wants to fork.
     *
     * The operation delegates the request to reserveFreeNodeForServer() or
     * reserveFreeNodeForClient() depending on the rank of the process. If the
     * reservation is triggered by the node pool process server,
     * replyToMessages() is called before the reservation is done. This is for
     * fairness: Otherwise the node pool process would fork pretty often, before
     * a remote reservation is taken into account. This is because remote
     * queries are usually only answered if the node pool is waiting for any
     * message.
     */
    int reserveFreeNode();

    /**
     * Wait for nodes in the cluster.
     *
     * Well, just wait until all the nodes in the cluster announce that they are
     * idle, i.e. are waiting for something to do. A node is set idle in the
     * central administration as soon as the node pool receives a job request
     * message. If it receives a job request message before the node has
     * registered, the node pool waits for the registration before it continues.
     * Consequently, this operation acts as a global barrier if it is called in
     * the beginning - a barrier that waits for two messages from each node.
     */
    void waitForAllNodesToBecomeIdle();

    /**
     * Master of This Worker Process
     *
     * For non-idle processes this operation returns the number of the process
     * that awakened them. If you want to find out the global master of the
     * cluster, you have to call Node::getMasterProcessRank().
     */
    int getMasterRank() const;

    /**
     * Handles all the message sent to the server.
     */
    virtual void receiveDanglingMessages();

    /**
     * Tell the pool to use a different strategy. If you pass such an object,
     * the node pool takes care to destroy it in the end (or if you reset it
     * later on). Please ensure that the node pool is shut down before you
     * reset the strategy, i.e. that you've called terminate and the
     * corresponding wait. After you've set a new strategy, call restart().
     */
    void setStrategy(NodePoolStrategy* strategy);

    /**
     * Tag for Fork Message
     *
     * Each worker listens on this tag for a fork message.
     *
     * Originall, I wanted ot use the tag afterwards also to exchange the
     * join and fork data. This is not a good idea, as the forks and joins
     * are realised within a service that polls the mpi queue permanently.
     * Consequently, it might happen that the fork message overtakes the
     * fork and join data and this leads to a break down of the algorithm.
     *
     * @return Tag for fork messages and fork and join data.
     */
    int getTagForForkMessages() const;

    /**
     * Wake up node pool or register at the node pool. Depends on where the
     * restart is called. restart() acts as barrier, i.e. ensure that you call
     * it on all nodes at the same time.
     */
    void restart();

    /**
     * Tells the node pool to shut down. There is a slight difference between
     * terminate() and shutdown(). The counterpart of shutdown() is init() and
     * both operations are called once throughout the whole application
     * lifetime. In contrast, terminate() and restart() are called several
     * times.
     */
    void shutdown();

    /**
     * This operation returns the number of non-idle nodes including the node
     * pool rank itself.
     */
    int getNumberOfWorkingNodes() const;

    /**
     * Broadcast to non-idle nodes
     *
     * Take an MPI message and broadcast it to all non-idle nodes that are
     * registered at the node pool. The message has to be a DaStGen object. You
     * are allowed to call this function only on the global rank 0. In return,
     * the broadcasted message never is sent to rank 0, i.e. the global master.
     *
     * Works on the master node only.
     */
    template <class Message>
    void broadcastToWorkingNodes(Message& message, int tag);

    void answerAllJobRequestMessagesWithHandleLocalProblem( bool value );

    /**
     * Send a job request message answer value ExecuteGlobalTask to all idle nodes
     * and nodes which have sent a job request at the moment.
     */
    void sendExecuteGlobalTaskJobMessages();
};


#include "tarch/parallel/NodePool.cpph"

#endif
