// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_PARALLEL_NODE_H_
#define _TARCH_PARALLEL_NODE_H_

#ifdef Parallel
#include <mpi.h>

#include "tarch/logging/Log.h"


namespace tarch {
  namespace parallel {
    class Node;

    /**
     * Returns a string representation of the mpi status. For a detailed
     * description of the mpi status itself see the file mpi.h.
     */
#ifdef Parallel
    std::string MPIStatusToString( const MPI_Status& status );

    std::string MPIReturnValueToString( int result );
#endif
  }
}

/**
 * Parallel Node Representant
 *
 * Represents a program instance within a cluster. Thus, this class is a
 * singleton.
 *
 * The parallel concept is a client - server model (process 0 is the server),
 * where all active nodes act as servers deploying calculations on demand. So
 * the basic activities of a parallel node are
 *
 * - receive new root element to work on
 * - pass back space-tree
 * - perform additive cycle
 * - perform multiplicative cycle
 *
 * The two perform commands have to return several statistic records. Among them
 * are time needed, some residual characteristics and flags defining if the
 * domain has been refined. Furthermore the number of communication partners is
 * an interesting detail.
 *
 * In the near future, this class should become responsible for the error
 * handling. Right now, the error handle set is the fatal error handler. That
 * is the whole parallel application is shut down as soon as an mpi error
 * occurs.
 *
 * @author Tobias Weinzierl
 * @version $Revision: 1.51 $
 */
class tarch::parallel::Node {
  public:
    static const int DEADLOCK_EXIT_CODE = -2;
  private:
    /**
     * Logging device.
     */
    static tarch::logging::Log _log;

    /**
     * Is set true if init() is called.
     */
    static bool _initIsCalled;

    /**
     * Rank (id) of this process.
     */
    int _rank;

    /**
     * Number of processors available.
     */
    int _numberOfProcessors;
#ifdef Parallel
    /**
     * MPI Communicator this process belongs to.
     */
    MPI_Comm _communicator;
#endif
    /**
     * How long shall the application wait until it writes a time-out warning
     */
    clock_t _timeOutWarning;

    clock_t _deadlockTimeOut;

    /**
     * The standard constructor assignes the attributes default values and
     * checks whether the program is compiled using the -DParallel option.
     * If this is not the case, a warning is logged.
     */
    Node();

    /**
     * The copy constructor is private.
     */
    Node( const Node& node );

    /**
     * Receive any Message Pending in the MPI/Receive Buffers
     */
    void receiveDanglingMessagesFromReceiveBuffers();

  public:
    /**
     * Return a Free Tag
     *
     * Returns a free tag to be used for a new datatype. Each result is
     * delivered exactly once. The string argument is just for logging.
     *
     *
     * !!! Implementation details
     *
     * This operation should write something to the log devices. However, it
     * is static and the class' log devices are static, too. C++ has no
     * mechanism to define which static entitiy has to be instantiated first.
     * On some systems, it hence happened that someone called this static
     * function while the static log attribute has not been initialised yet.
     * Consequently, the operation uses its own (local) log variable instead
     * of the log variable of the class. This is an important workaround!.
     */
    static int reserveFreeTag(const std::string& fullQualifiedMessageName);

    /**
     * Logs the status of the process onto the log device.
     *
     */
    void logStatus() const;

    /**
     * This operation returns the singleton instance. Before using this
     * instance, one has to call the init() operation on the instance returned.
     *
     * @return The singleton instance
     */
    static Node& getInstance();

    /**
     * The standard destructor calls MPI_Finalize().
     */
    virtual ~Node();

    /**
     * This operation initializes the MPI environment and the program instance.
     * Not that the argv and argc parameters are both in and out parameters.
     * Before you pass them to the operation, they are set by the mpi
     * environment. Afterwards the original parameters from the user are stored
     * within them.
     *
     * !! Implementation details
     *
     * init never uses the log device to report any errors as the log device
     * usually in turn uses Node's getters. Furthermore, the _initIsCalled flag
     * has thus to be set before the log state operation is invoked.
     *
     * @return true if initialisation has been successful
     */
    bool init(int* argc, char*** argv);

    /**
     * Shuts down the application. Should be the last operation called by the
     * overall application.
     *
     * !!! Rationale
     *
     * Originally, I put the shutdown operation into the destructor of Node.
     * The MPI environment consequently was shut down as soon as the operating
     * system terminates the application. However, Scalasca complained on the
     * BlueGene/P that the
     * destruction happened after the return statement of the main method. To
     * make Peano work with Scalasca, I hence moved the MPI shutdown into a
     * method called explicitely before the final return statement.
     *
     * This seems to be a CLX effect, i.e. the Intel and the GNU compilers
     * worked fine with Scalasca. Hence, I assume that Intel and GNU
     * executables destroy all static objects (singletons) before they return
     * from the main functions. CLX destroys the static objects after the
     * return statement and thus makes Scalasca's instrumentation report an
     * error.
     */
    void shutdown();

    /**
     * @return Rank of this node
     */
    int getRank() const;

    /**
     * @return 0
     */
    static int getGlobalMasterRank();
#ifdef Parallel
    /**
     * @return Communicator of Peano
     */
    MPI_Comm getCommunicator() const;
#endif
    /**
     * @return Number of Nodes Available
     */
    int getNumberOfNodes() const;

    /**
     * @return Is this node the global master process, i.e. does its rank equal getMasterProcessRank()?
     */
    bool isGlobalMaster() const;

    /**
     * Triggers a time out and shuts down the cluster:
     *
     * The implementation does not use MPI_Abort, since it seems that this
     * operation requires all nodes running. Instead of,
     * getDeadlockWarningTimeStamp() uses the system exit function passing it
     * DEADLOCK_EXIT_CODE as exit code.
     *
     * The operation should be called only if the deadlock time-out is switched
     * on ( isTimeOutDeadlockEnabled() ) and the deadlock time-out has expired.
     * Use getDeadlockWarningTimeStamp() and the system operation clock() to
     * check the second requirement.
     *
     * @param className  Name of the class that triggers the deadlock shutdown.
     * @param methodName Name of the method that triggers the deadlock shutdown.
     * @param communicationPartnerRank Rank of the node the operation that
     *          should have sent a message but did not.
     */
    void triggerDeadlockTimeOut(
      const std::string&  className,
      const std::string&  methodName,
      int                 communicationPartnerRank,
      int                 tag,
      int                 numberOfExpectedMessages,
      const std::string&  comment = ""
    );

    void writeTimeOutWarning(
      const std::string&  className,
      const std::string&  methodName,
      int                 communicationPartnerRank,
      int                 tag,
      int                 numberOfExpectedMessages
    );

    /**
     * Old signature calls new function and gives warning that old signature is used.
     */
    void triggerDeadlockTimeOut(
      const std::string&  className,
      const std::string&  methodName,
      int                 communicationPartnerRank,
      const std::string&  comment = ""
    );

    /**
     * Old signature calls new function and gives warning that old signature is used.
     */
    void writeTimeOutWarning(
      const std::string&  className,
      const std::string&  methodName,
      int                 communicationPartnerRank
    );

    void plotMessageQueues();

    /**
     * Ensure that there are no messages anymore from the specified rank.
     */
    void ensureThatMessageQueuesAreEmpty( int fromRank, int tag );

    /**
     * If you want to make any receive operation to be able to recognize a
     * time-out, you have to insert the following code snipped before your
     * receive call.
     *
     * \code
     *   MPI_Status  status;
     *   int         flag = 0;
     *   clock_t     timeOutWarning   = Node::getInstance().getDeadlockWarningTimeStamp();
     *   clock_t     timeOutShutdown  = Node::getInstance().getDeadlockTimeOutTimeStamp();
     *   bool        triggeredTimeoutWarning = false;
     *   MPI_Iprobe( source, CommunicationTag, Node::getInstance().getCommunicator(), &flag, &status );
     *
     *   while (!flag) {
     *     if ( Node::getInstance().isTimeOutWarningEnabled() && (clock()>timeOutWarning) && (!triggeredTimeoutWarning)) {
     *       Node::getInstance().writeTimeOutWarning( "parallel::SendReceiveBuffer", "receiveAllMessages()", source );
     *       triggeredTimeoutWarning = true;
     *     }
     *     if ( Node::getInstance().isTimeOutDeadlockEnabled() && (clock()>timeOutShutdown)) {
     *       Node::getInstance().triggerDeadlockTimeOut( "parallel::SendReceiveBuffer", "receiveAllMessages()", source );
     *     }
     *     MPI_Iprobe( source, CommunicationTag, Node::getInstance().getCommunicator(), &flag, &status );
     *   }
     * \endcode
     *
     * @return Time stamp when next warning should be written if no message has
     *         been received meanwhile.
     */
    clock_t getDeadlockWarningTimeStamp() const;

    /**
     * @return Time stamp when next application should terminate because of a
     *         time out if no message has been received meanwhile.
     */
    clock_t getDeadlockTimeOutTimeStamp() const;

    bool isTimeOutDeadlockEnabled() const;
    bool isTimeOutWarningEnabled() const;

    bool isInitialised() const;

    void setTimeOutWarning( const clock_t & value );
    void setDeadlockTimeOut( const clock_t & value );
#ifdef Parallel
    void setCommunicator( MPI_Comm communicator );
#endif

    /**
     * Receive Dangling MPI Messages
     *
     * This operation tells all the services that they shall poll the MPI queues
     * for additional messages. Thus, it is questionable why such an operation
     * is assigned to this class. It only wraps the services. However, all these
     * 'hey poll the MPI queues' is always triggered by blocking sends and
     * receives, i.e. it is tied to the parallelisation. And, thus, it is part
     * of the node singleton.
     */
    void receiveDanglingMessages();


};
#endif
#endif
