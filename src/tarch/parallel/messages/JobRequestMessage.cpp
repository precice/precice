#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/parallel/messages/JobRequestMessage.h"

tarch::parallel::messages::JobRequestMessage::PersistentRecords::PersistentRecords() {
   
}


tarch::parallel::messages::JobRequestMessage::PersistentRecords::PersistentRecords(const bool& dummy):
_dummy(dummy) {
   
}


 bool tarch::parallel::messages::JobRequestMessage::PersistentRecords::getDummy() const  {
   return _dummy;
}



 void tarch::parallel::messages::JobRequestMessage::PersistentRecords::setDummy(const bool& dummy)  {
   _dummy = dummy;
}


tarch::parallel::messages::JobRequestMessage::JobRequestMessage() {
   
}


tarch::parallel::messages::JobRequestMessage::JobRequestMessage(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._dummy) {
   
}


tarch::parallel::messages::JobRequestMessage::JobRequestMessage(const bool& dummy):
_persistentRecords(dummy) {
   
}


tarch::parallel::messages::JobRequestMessage::~JobRequestMessage() { }


 bool tarch::parallel::messages::JobRequestMessage::getDummy() const  {
   return _persistentRecords._dummy;
}



 void tarch::parallel::messages::JobRequestMessage::setDummy(const bool& dummy)  {
   _persistentRecords._dummy = dummy;
}




std::string tarch::parallel::messages::JobRequestMessage::toString() const {
   std::ostringstream stringstr;
   toString(stringstr);
   return stringstr.str();
}

void tarch::parallel::messages::JobRequestMessage::toString (std::ostream& out) const {
   out << "("; 
   out << "dummy:" << getDummy();
   out <<  ")";
}


tarch::parallel::messages::JobRequestMessage::PersistentRecords tarch::parallel::messages::JobRequestMessage::getPersistentRecords() const {
   return _persistentRecords;
}

tarch::parallel::messages::JobRequestMessagePacked tarch::parallel::messages::JobRequestMessage::convert() const{
   return JobRequestMessagePacked(
      getDummy()
   );
}

#ifdef Parallel
   tarch::logging::Log tarch::parallel::messages::JobRequestMessage::_log( "tarch::parallel::messages::JobRequestMessage" );
   
   MPI_Datatype tarch::parallel::messages::JobRequestMessage::Datatype = 0;
   MPI_Datatype tarch::parallel::messages::JobRequestMessage::FullDatatype = 0;
   
   
   void tarch::parallel::messages::JobRequestMessage::initDatatype() {
      {
         JobRequestMessage dummyJobRequestMessage[2];
         
         const int Attributes = 2;
         MPI_Datatype subtypes[Attributes] = {
            MPI_CHAR,		 //dummy
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //dummy
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyJobRequestMessage[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyJobRequestMessage[0]._persistentRecords._dummy))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyJobRequestMessage[1]._persistentRecords._dummy))), 		&disp[1] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &JobRequestMessage::Datatype );
         MPI_Type_commit( &JobRequestMessage::Datatype );
         
      }
      {
         JobRequestMessage dummyJobRequestMessage[2];
         
         const int Attributes = 2;
         MPI_Datatype subtypes[Attributes] = {
            MPI_CHAR,		 //dummy
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //dummy
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyJobRequestMessage[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyJobRequestMessage[0]._persistentRecords._dummy))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyJobRequestMessage[1]._persistentRecords._dummy))), 		&disp[1] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &JobRequestMessage::FullDatatype );
         MPI_Type_commit( &JobRequestMessage::FullDatatype );
         
      }
      
   }
   
   
   void tarch::parallel::messages::JobRequestMessage::shutdownDatatype() {
      MPI_Type_free( &JobRequestMessage::Datatype );
      MPI_Type_free( &JobRequestMessage::FullDatatype );
      
   }
   
   void tarch::parallel::messages::JobRequestMessage::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
         msg << "was not able to send message tarch::parallel::messages::JobRequestMessage "
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
            msg << "testing for finished send task for tarch::parallel::messages::JobRequestMessage "
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
            "tarch::parallel::messages::JobRequestMessage",
            "send(int)", destination,tag,1
            );
            triggeredTimeoutWarning = true;
         }
         if (
            tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
            (clock()>timeOutShutdown)
         ) {
            tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "tarch::parallel::messages::JobRequestMessage",
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
   
   
   
   void tarch::parallel::messages::JobRequestMessage::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
         msg << "failed to start to receive tarch::parallel::messages::JobRequestMessage from node "
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
            msg << "testing for finished receive task for tarch::parallel::messages::JobRequestMessage failed: "
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
            "tarch::parallel::messages::JobRequestMessage",
            "receive(int)", source,tag,1
            );
            triggeredTimeoutWarning = true;
         }
         if (
            tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
            (clock()>timeOutShutdown)
         ) {
            tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "tarch::parallel::messages::JobRequestMessage",
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
   
   
   
   bool tarch::parallel::messages::JobRequestMessage::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
   
   int tarch::parallel::messages::JobRequestMessage::getSenderRank() const {
      assertion( _senderRank!=-1 );
      return _senderRank;
      
   }
#endif


tarch::parallel::messages::JobRequestMessagePacked::PersistentRecords::PersistentRecords() {
   assertion((1 < (8 * sizeof(short int))));
   
}


tarch::parallel::messages::JobRequestMessagePacked::PersistentRecords::PersistentRecords(const bool& dummy) {
   setDummy(dummy);
   assertion((1 < (8 * sizeof(short int))));
   
}


 bool tarch::parallel::messages::JobRequestMessagePacked::PersistentRecords::getDummy() const  {
   short int mask = 1 << (0);
   short int tmp = static_cast<short int>(_packedRecords0 & mask);
   return (tmp != 0);
}



 void tarch::parallel::messages::JobRequestMessagePacked::PersistentRecords::setDummy(const bool& dummy)  {
   short int mask = 1 << (0);
   _packedRecords0 = static_cast<short int>( dummy ? (_packedRecords0 | mask) : (_packedRecords0 & ~mask));
}


tarch::parallel::messages::JobRequestMessagePacked::JobRequestMessagePacked() {
   assertion((1 < (8 * sizeof(short int))));
   
}


tarch::parallel::messages::JobRequestMessagePacked::JobRequestMessagePacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords.getDummy()) {
   assertion((1 < (8 * sizeof(short int))));
   
}


tarch::parallel::messages::JobRequestMessagePacked::JobRequestMessagePacked(const bool& dummy):
_persistentRecords(dummy) {
   assertion((1 < (8 * sizeof(short int))));
   
}


tarch::parallel::messages::JobRequestMessagePacked::~JobRequestMessagePacked() { }


 bool tarch::parallel::messages::JobRequestMessagePacked::getDummy() const  {
   short int mask = 1 << (0);
   short int tmp = static_cast<short int>(_persistentRecords._packedRecords0 & mask);
   return (tmp != 0);
}



 void tarch::parallel::messages::JobRequestMessagePacked::setDummy(const bool& dummy)  {
   short int mask = 1 << (0);
   _persistentRecords._packedRecords0 = static_cast<short int>( dummy ? (_persistentRecords._packedRecords0 | mask) : (_persistentRecords._packedRecords0 & ~mask));
}




std::string tarch::parallel::messages::JobRequestMessagePacked::toString() const {
   std::ostringstream stringstr;
   toString(stringstr);
   return stringstr.str();
}

void tarch::parallel::messages::JobRequestMessagePacked::toString (std::ostream& out) const {
   out << "("; 
   out << "dummy:" << getDummy();
   out <<  ")";
}


tarch::parallel::messages::JobRequestMessagePacked::PersistentRecords tarch::parallel::messages::JobRequestMessagePacked::getPersistentRecords() const {
   return _persistentRecords;
}

tarch::parallel::messages::JobRequestMessage tarch::parallel::messages::JobRequestMessagePacked::convert() const{
   return JobRequestMessage(
      getDummy()
   );
}

#ifdef Parallel
   tarch::logging::Log tarch::parallel::messages::JobRequestMessagePacked::_log( "tarch::parallel::messages::JobRequestMessagePacked" );
   
   MPI_Datatype tarch::parallel::messages::JobRequestMessagePacked::Datatype = 0;
   MPI_Datatype tarch::parallel::messages::JobRequestMessagePacked::FullDatatype = 0;
   
   
   void tarch::parallel::messages::JobRequestMessagePacked::initDatatype() {
      {
         JobRequestMessagePacked dummyJobRequestMessagePacked[2];
         
         const int Attributes = 2;
         MPI_Datatype subtypes[Attributes] = {
            MPI_SHORT,		 //_packedRecords0
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //_packedRecords0
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyJobRequestMessagePacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyJobRequestMessagePacked[0]._persistentRecords._packedRecords0))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyJobRequestMessagePacked[1]._persistentRecords._packedRecords0))), 		&disp[1] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &JobRequestMessagePacked::Datatype );
         MPI_Type_commit( &JobRequestMessagePacked::Datatype );
         
      }
      {
         JobRequestMessagePacked dummyJobRequestMessagePacked[2];
         
         const int Attributes = 2;
         MPI_Datatype subtypes[Attributes] = {
            MPI_SHORT,		 //_packedRecords0
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //_packedRecords0
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyJobRequestMessagePacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyJobRequestMessagePacked[0]._persistentRecords._packedRecords0))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyJobRequestMessagePacked[1]._persistentRecords._packedRecords0))), 		&disp[1] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &JobRequestMessagePacked::FullDatatype );
         MPI_Type_commit( &JobRequestMessagePacked::FullDatatype );
         
      }
      
   }
   
   
   void tarch::parallel::messages::JobRequestMessagePacked::shutdownDatatype() {
      MPI_Type_free( &JobRequestMessagePacked::Datatype );
      MPI_Type_free( &JobRequestMessagePacked::FullDatatype );
      
   }
   
   void tarch::parallel::messages::JobRequestMessagePacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
         msg << "was not able to send message tarch::parallel::messages::JobRequestMessagePacked "
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
            msg << "testing for finished send task for tarch::parallel::messages::JobRequestMessagePacked "
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
            "tarch::parallel::messages::JobRequestMessagePacked",
            "send(int)", destination,tag,1
            );
            triggeredTimeoutWarning = true;
         }
         if (
            tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
            (clock()>timeOutShutdown)
         ) {
            tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "tarch::parallel::messages::JobRequestMessagePacked",
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
   
   
   
   void tarch::parallel::messages::JobRequestMessagePacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
         msg << "failed to start to receive tarch::parallel::messages::JobRequestMessagePacked from node "
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
            msg << "testing for finished receive task for tarch::parallel::messages::JobRequestMessagePacked failed: "
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
            "tarch::parallel::messages::JobRequestMessagePacked",
            "receive(int)", source,tag,1
            );
            triggeredTimeoutWarning = true;
         }
         if (
            tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
            (clock()>timeOutShutdown)
         ) {
            tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "tarch::parallel::messages::JobRequestMessagePacked",
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
   
   
   
   bool tarch::parallel::messages::JobRequestMessagePacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
   
   int tarch::parallel::messages::JobRequestMessagePacked::getSenderRank() const {
      assertion( _senderRank!=-1 );
      return _senderRank;
      
   }
#endif



