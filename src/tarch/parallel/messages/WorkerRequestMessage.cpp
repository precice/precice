#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/parallel/messages/WorkerRequestMessage.h"

tarch::parallel::messages::WorkerRequestMessage::PersistentRecords::PersistentRecords() {
   
}


tarch::parallel::messages::WorkerRequestMessage::PersistentRecords::PersistentRecords(const int& tmp):
_tmp(tmp) {
   
}


 int tarch::parallel::messages::WorkerRequestMessage::PersistentRecords::getTmp() const  {
   return _tmp;
}



 void tarch::parallel::messages::WorkerRequestMessage::PersistentRecords::setTmp(const int& tmp)  {
   _tmp = tmp;
}


tarch::parallel::messages::WorkerRequestMessage::WorkerRequestMessage() {
   
}


tarch::parallel::messages::WorkerRequestMessage::WorkerRequestMessage(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._tmp) {
   
}


tarch::parallel::messages::WorkerRequestMessage::WorkerRequestMessage(const int& tmp):
_persistentRecords(tmp) {
   
}


tarch::parallel::messages::WorkerRequestMessage::~WorkerRequestMessage() { }


 int tarch::parallel::messages::WorkerRequestMessage::getTmp() const  {
   return _persistentRecords._tmp;
}



 void tarch::parallel::messages::WorkerRequestMessage::setTmp(const int& tmp)  {
   _persistentRecords._tmp = tmp;
}




std::string tarch::parallel::messages::WorkerRequestMessage::toString() const {
   std::ostringstream stringstr;
   toString(stringstr);
   return stringstr.str();
}

void tarch::parallel::messages::WorkerRequestMessage::toString (std::ostream& out) const {
   out << "("; 
   out << "tmp:" << getTmp();
   out <<  ")";
}


tarch::parallel::messages::WorkerRequestMessage::PersistentRecords tarch::parallel::messages::WorkerRequestMessage::getPersistentRecords() const {
   return _persistentRecords;
}

tarch::parallel::messages::WorkerRequestMessagePacked tarch::parallel::messages::WorkerRequestMessage::convert() const{
   return WorkerRequestMessagePacked(
      getTmp()
   );
}

#ifdef Parallel
   tarch::logging::Log tarch::parallel::messages::WorkerRequestMessage::_log( "tarch::parallel::messages::WorkerRequestMessage" );
   
   MPI_Datatype tarch::parallel::messages::WorkerRequestMessage::Datatype = 0;
   MPI_Datatype tarch::parallel::messages::WorkerRequestMessage::FullDatatype = 0;
   
   
   void tarch::parallel::messages::WorkerRequestMessage::initDatatype() {
      {
         WorkerRequestMessage dummyWorkerRequestMessage[2];
         
         const int Attributes = 2;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //tmp
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //tmp
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyWorkerRequestMessage[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyWorkerRequestMessage[0]._persistentRecords._tmp))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyWorkerRequestMessage[1]._persistentRecords._tmp))), 		&disp[1] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &WorkerRequestMessage::Datatype );
         MPI_Type_commit( &WorkerRequestMessage::Datatype );
         
      }
      {
         WorkerRequestMessage dummyWorkerRequestMessage[2];
         
         const int Attributes = 2;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //tmp
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //tmp
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyWorkerRequestMessage[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyWorkerRequestMessage[0]._persistentRecords._tmp))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyWorkerRequestMessage[1]._persistentRecords._tmp))), 		&disp[1] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &WorkerRequestMessage::FullDatatype );
         MPI_Type_commit( &WorkerRequestMessage::FullDatatype );
         
      }
      
   }
   
   
   void tarch::parallel::messages::WorkerRequestMessage::shutdownDatatype() {
      MPI_Type_free( &WorkerRequestMessage::Datatype );
      MPI_Type_free( &WorkerRequestMessage::FullDatatype );
      
   }
   
   void tarch::parallel::messages::WorkerRequestMessage::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
      MPI_Request* sendRequestHandle = new MPI_Request();
      MPI_Status   status;
      int          flag = 0;
      int          result;
      
      clock_t      timeOutWarning   = -1;
      clock_t      timeOutShutdown  = -1;
      bool         triggeredTimeoutWarning = false;
      
      #ifdef Asserts
      _senderRank = -1;
      #endif
      
      if (exchangeOnlyAttributesMarkedWithParallelise) {
         result = MPI_Isend(
            this, 1, Datatype, destination,
            tag, tarch::parallel::Node::getInstance().getCommunicator(),
            sendRequestHandle
         );
         
      }
      else {
         result = MPI_Isend(
            this, 1, FullDatatype, destination,
            tag, tarch::parallel::Node::getInstance().getCommunicator(),
            sendRequestHandle
         );
         
      }
      if  (result!=MPI_SUCCESS) {
         std::ostringstream msg;
         msg << "was not able to send message tarch::parallel::messages::WorkerRequestMessage "
         << toString()
         << " to node " << destination
         << ": " << tarch::parallel::MPIReturnValueToString(result);
         _log.error( "send(int)",msg.str() );
      }
      result = MPI_Test( sendRequestHandle, &flag, &status );
      while (!flag) {
         if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
         if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
         result = MPI_Test( sendRequestHandle, &flag, &status );
         if (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "testing for finished send task for tarch::parallel::messages::WorkerRequestMessage "
            << toString()
            << " sent to node " << destination
            << " failed: " << tarch::parallel::MPIReturnValueToString(result);
            _log.error("send(int)", msg.str() );
         }
         
         // deadlock aspect
         if (
            tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
            (clock()>timeOutWarning) &&
            (!triggeredTimeoutWarning)
         ) {
            tarch::parallel::Node::getInstance().writeTimeOutWarning(
            "tarch::parallel::messages::WorkerRequestMessage",
            "send(int)", destination,tag,1
            );
            triggeredTimeoutWarning = true;
         }
         if (
            tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
            (clock()>timeOutShutdown)
         ) {
            tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "tarch::parallel::messages::WorkerRequestMessage",
            "send(int)", destination,tag,1
            );
         }
         tarch::parallel::Node::getInstance().receiveDanglingMessages();
      }
      
      delete sendRequestHandle;
      #ifdef Debug
      _log.debug("send(int,int)", "sent " + toString() );
      #endif
      
   }
   
   
   
   void tarch::parallel::messages::WorkerRequestMessage::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
      MPI_Request* sendRequestHandle = new MPI_Request();
      MPI_Status   status;
      int          flag = 0;
      int          result;
      
      clock_t      timeOutWarning   = -1;
      clock_t      timeOutShutdown  = -1;
      bool         triggeredTimeoutWarning = false;
      
      if (exchangeOnlyAttributesMarkedWithParallelise) {
         result = MPI_Irecv(
            this, 1, Datatype, source, tag,
            tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
         );
         
      }
      else {
         result = MPI_Irecv(
            this, 1, FullDatatype, source, tag,
            tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
         );
         
      }
      if ( result != MPI_SUCCESS ) {
         std::ostringstream msg;
         msg << "failed to start to receive tarch::parallel::messages::WorkerRequestMessage from node "
         << source << ": " << tarch::parallel::MPIReturnValueToString(result);
         _log.error( "receive(int)", msg.str() );
      }
      
      result = MPI_Test( sendRequestHandle, &flag, &status );
      while (!flag) {
         if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
         if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
         result = MPI_Test( sendRequestHandle, &flag, &status );
         if (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "testing for finished receive task for tarch::parallel::messages::WorkerRequestMessage failed: "
            << tarch::parallel::MPIReturnValueToString(result);
            _log.error("receive(int)", msg.str() );
         }
         
         // deadlock aspect
         if (
            tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
            (clock()>timeOutWarning) &&
            (!triggeredTimeoutWarning)
         ) {
            tarch::parallel::Node::getInstance().writeTimeOutWarning(
            "tarch::parallel::messages::WorkerRequestMessage",
            "receive(int)", source,tag,1
            );
            triggeredTimeoutWarning = true;
         }
         if (
            tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
            (clock()>timeOutShutdown)
         ) {
            tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "tarch::parallel::messages::WorkerRequestMessage",
            "receive(int)", source,tag,1
            );
         }
         tarch::parallel::Node::getInstance().receiveDanglingMessages();
      }
      
      delete sendRequestHandle;
      
      _senderRank = status.MPI_SOURCE;
      #ifdef Debug
      _log.debug("receive(int,int)", "received " + toString() ); 
      #endif
      
   }
   
   
   
   bool tarch::parallel::messages::WorkerRequestMessage::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
      MPI_Status status;
      int  flag        = 0;
      MPI_Iprobe(
         MPI_ANY_SOURCE, tag,
         tarch::parallel::Node::getInstance().getCommunicator(), &flag, &status
      );
      if (flag) {
         int  messageCounter;
         if (exchangeOnlyAttributesMarkedWithParallelise) {
            MPI_Get_count(&status, Datatype, &messageCounter);
         }
         else {
            MPI_Get_count(&status, FullDatatype, &messageCounter);
         }
         return messageCounter > 0;
      }
      else return false;
      
   }
   
   int tarch::parallel::messages::WorkerRequestMessage::getSenderRank() const {
      assertion( _senderRank!=-1 );
      return _senderRank;
      
   }
#endif


tarch::parallel::messages::WorkerRequestMessagePacked::PersistentRecords::PersistentRecords() {
   
}


tarch::parallel::messages::WorkerRequestMessagePacked::PersistentRecords::PersistentRecords(const int& tmp):
_tmp(tmp) {
   
}


 int tarch::parallel::messages::WorkerRequestMessagePacked::PersistentRecords::getTmp() const  {
   return _tmp;
}



 void tarch::parallel::messages::WorkerRequestMessagePacked::PersistentRecords::setTmp(const int& tmp)  {
   _tmp = tmp;
}


tarch::parallel::messages::WorkerRequestMessagePacked::WorkerRequestMessagePacked() {
   
}


tarch::parallel::messages::WorkerRequestMessagePacked::WorkerRequestMessagePacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._tmp) {
   
}


tarch::parallel::messages::WorkerRequestMessagePacked::WorkerRequestMessagePacked(const int& tmp):
_persistentRecords(tmp) {
   
}


tarch::parallel::messages::WorkerRequestMessagePacked::~WorkerRequestMessagePacked() { }


 int tarch::parallel::messages::WorkerRequestMessagePacked::getTmp() const  {
   return _persistentRecords._tmp;
}



 void tarch::parallel::messages::WorkerRequestMessagePacked::setTmp(const int& tmp)  {
   _persistentRecords._tmp = tmp;
}




std::string tarch::parallel::messages::WorkerRequestMessagePacked::toString() const {
   std::ostringstream stringstr;
   toString(stringstr);
   return stringstr.str();
}

void tarch::parallel::messages::WorkerRequestMessagePacked::toString (std::ostream& out) const {
   out << "("; 
   out << "tmp:" << getTmp();
   out <<  ")";
}


tarch::parallel::messages::WorkerRequestMessagePacked::PersistentRecords tarch::parallel::messages::WorkerRequestMessagePacked::getPersistentRecords() const {
   return _persistentRecords;
}

tarch::parallel::messages::WorkerRequestMessage tarch::parallel::messages::WorkerRequestMessagePacked::convert() const{
   return WorkerRequestMessage(
      getTmp()
   );
}

#ifdef Parallel
   tarch::logging::Log tarch::parallel::messages::WorkerRequestMessagePacked::_log( "tarch::parallel::messages::WorkerRequestMessagePacked" );
   
   MPI_Datatype tarch::parallel::messages::WorkerRequestMessagePacked::Datatype = 0;
   MPI_Datatype tarch::parallel::messages::WorkerRequestMessagePacked::FullDatatype = 0;
   
   
   void tarch::parallel::messages::WorkerRequestMessagePacked::initDatatype() {
      {
         WorkerRequestMessagePacked dummyWorkerRequestMessagePacked[2];
         
         const int Attributes = 2;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //tmp
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //tmp
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyWorkerRequestMessagePacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyWorkerRequestMessagePacked[0]._persistentRecords._tmp))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyWorkerRequestMessagePacked[1]._persistentRecords._tmp))), 		&disp[1] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &WorkerRequestMessagePacked::Datatype );
         MPI_Type_commit( &WorkerRequestMessagePacked::Datatype );
         
      }
      {
         WorkerRequestMessagePacked dummyWorkerRequestMessagePacked[2];
         
         const int Attributes = 2;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //tmp
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //tmp
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyWorkerRequestMessagePacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyWorkerRequestMessagePacked[0]._persistentRecords._tmp))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyWorkerRequestMessagePacked[1]._persistentRecords._tmp))), 		&disp[1] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &WorkerRequestMessagePacked::FullDatatype );
         MPI_Type_commit( &WorkerRequestMessagePacked::FullDatatype );
         
      }
      
   }
   
   
   void tarch::parallel::messages::WorkerRequestMessagePacked::shutdownDatatype() {
      MPI_Type_free( &WorkerRequestMessagePacked::Datatype );
      MPI_Type_free( &WorkerRequestMessagePacked::FullDatatype );
      
   }
   
   void tarch::parallel::messages::WorkerRequestMessagePacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
      MPI_Request* sendRequestHandle = new MPI_Request();
      MPI_Status   status;
      int          flag = 0;
      int          result;
      
      clock_t      timeOutWarning   = -1;
      clock_t      timeOutShutdown  = -1;
      bool         triggeredTimeoutWarning = false;
      
      #ifdef Asserts
      _senderRank = -1;
      #endif
      
      if (exchangeOnlyAttributesMarkedWithParallelise) {
         result = MPI_Isend(
            this, 1, Datatype, destination,
            tag, tarch::parallel::Node::getInstance().getCommunicator(),
            sendRequestHandle
         );
         
      }
      else {
         result = MPI_Isend(
            this, 1, FullDatatype, destination,
            tag, tarch::parallel::Node::getInstance().getCommunicator(),
            sendRequestHandle
         );
         
      }
      if  (result!=MPI_SUCCESS) {
         std::ostringstream msg;
         msg << "was not able to send message tarch::parallel::messages::WorkerRequestMessagePacked "
         << toString()
         << " to node " << destination
         << ": " << tarch::parallel::MPIReturnValueToString(result);
         _log.error( "send(int)",msg.str() );
      }
      result = MPI_Test( sendRequestHandle, &flag, &status );
      while (!flag) {
         if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
         if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
         result = MPI_Test( sendRequestHandle, &flag, &status );
         if (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "testing for finished send task for tarch::parallel::messages::WorkerRequestMessagePacked "
            << toString()
            << " sent to node " << destination
            << " failed: " << tarch::parallel::MPIReturnValueToString(result);
            _log.error("send(int)", msg.str() );
         }
         
         // deadlock aspect
         if (
            tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
            (clock()>timeOutWarning) &&
            (!triggeredTimeoutWarning)
         ) {
            tarch::parallel::Node::getInstance().writeTimeOutWarning(
            "tarch::parallel::messages::WorkerRequestMessagePacked",
            "send(int)", destination,tag,1
            );
            triggeredTimeoutWarning = true;
         }
         if (
            tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
            (clock()>timeOutShutdown)
         ) {
            tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "tarch::parallel::messages::WorkerRequestMessagePacked",
            "send(int)", destination,tag,1
            );
         }
         tarch::parallel::Node::getInstance().receiveDanglingMessages();
      }
      
      delete sendRequestHandle;
      #ifdef Debug
      _log.debug("send(int,int)", "sent " + toString() );
      #endif
      
   }
   
   
   
   void tarch::parallel::messages::WorkerRequestMessagePacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
      MPI_Request* sendRequestHandle = new MPI_Request();
      MPI_Status   status;
      int          flag = 0;
      int          result;
      
      clock_t      timeOutWarning   = -1;
      clock_t      timeOutShutdown  = -1;
      bool         triggeredTimeoutWarning = false;
      
      if (exchangeOnlyAttributesMarkedWithParallelise) {
         result = MPI_Irecv(
            this, 1, Datatype, source, tag,
            tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
         );
         
      }
      else {
         result = MPI_Irecv(
            this, 1, FullDatatype, source, tag,
            tarch::parallel::Node::getInstance().getCommunicator(), sendRequestHandle
         );
         
      }
      if ( result != MPI_SUCCESS ) {
         std::ostringstream msg;
         msg << "failed to start to receive tarch::parallel::messages::WorkerRequestMessagePacked from node "
         << source << ": " << tarch::parallel::MPIReturnValueToString(result);
         _log.error( "receive(int)", msg.str() );
      }
      
      result = MPI_Test( sendRequestHandle, &flag, &status );
      while (!flag) {
         if (timeOutWarning==-1)   timeOutWarning   = tarch::parallel::Node::getInstance().getDeadlockWarningTimeStamp();
         if (timeOutShutdown==-1)  timeOutShutdown  = tarch::parallel::Node::getInstance().getDeadlockTimeOutTimeStamp();
         result = MPI_Test( sendRequestHandle, &flag, &status );
         if (result!=MPI_SUCCESS) {
            std::ostringstream msg;
            msg << "testing for finished receive task for tarch::parallel::messages::WorkerRequestMessagePacked failed: "
            << tarch::parallel::MPIReturnValueToString(result);
            _log.error("receive(int)", msg.str() );
         }
         
         // deadlock aspect
         if (
            tarch::parallel::Node::getInstance().isTimeOutWarningEnabled() &&
            (clock()>timeOutWarning) &&
            (!triggeredTimeoutWarning)
         ) {
            tarch::parallel::Node::getInstance().writeTimeOutWarning(
            "tarch::parallel::messages::WorkerRequestMessagePacked",
            "receive(int)", source,tag,1
            );
            triggeredTimeoutWarning = true;
         }
         if (
            tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
            (clock()>timeOutShutdown)
         ) {
            tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "tarch::parallel::messages::WorkerRequestMessagePacked",
            "receive(int)", source,tag,1
            );
         }
         tarch::parallel::Node::getInstance().receiveDanglingMessages();
      }
      
      delete sendRequestHandle;
      
      _senderRank = status.MPI_SOURCE;
      #ifdef Debug
      _log.debug("receive(int,int)", "received " + toString() ); 
      #endif
      
   }
   
   
   
   bool tarch::parallel::messages::WorkerRequestMessagePacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
      MPI_Status status;
      int  flag        = 0;
      MPI_Iprobe(
         MPI_ANY_SOURCE, tag,
         tarch::parallel::Node::getInstance().getCommunicator(), &flag, &status
      );
      if (flag) {
         int  messageCounter;
         if (exchangeOnlyAttributesMarkedWithParallelise) {
            MPI_Get_count(&status, Datatype, &messageCounter);
         }
         else {
            MPI_Get_count(&status, FullDatatype, &messageCounter);
         }
         return messageCounter > 0;
      }
      else return false;
      
   }
   
   int tarch::parallel::messages::WorkerRequestMessagePacked::getSenderRank() const {
      assertion( _senderRank!=-1 );
      return _senderRank;
      
   }
#endif



