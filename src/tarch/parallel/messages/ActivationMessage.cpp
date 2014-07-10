#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/parallel/messages/ActivationMessage.h"

tarch::parallel::messages::ActivationMessage::PersistentRecords::PersistentRecords() {
   
}


tarch::parallel::messages::ActivationMessage::PersistentRecords::PersistentRecords(const int& newMaster):
_newMaster(newMaster) {
   
}


 int tarch::parallel::messages::ActivationMessage::PersistentRecords::getNewMaster() const  {
   return _newMaster;
}



 void tarch::parallel::messages::ActivationMessage::PersistentRecords::setNewMaster(const int& newMaster)  {
   _newMaster = newMaster;
}


tarch::parallel::messages::ActivationMessage::ActivationMessage() {
   
}


tarch::parallel::messages::ActivationMessage::ActivationMessage(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._newMaster) {
   
}


tarch::parallel::messages::ActivationMessage::ActivationMessage(const int& newMaster):
_persistentRecords(newMaster) {
   
}


tarch::parallel::messages::ActivationMessage::~ActivationMessage() { }


 int tarch::parallel::messages::ActivationMessage::getNewMaster() const  {
   return _persistentRecords._newMaster;
}



 void tarch::parallel::messages::ActivationMessage::setNewMaster(const int& newMaster)  {
   _persistentRecords._newMaster = newMaster;
}




std::string tarch::parallel::messages::ActivationMessage::toString() const {
   std::ostringstream stringstr;
   toString(stringstr);
   return stringstr.str();
}

void tarch::parallel::messages::ActivationMessage::toString (std::ostream& out) const {
   out << "("; 
   out << "newMaster:" << getNewMaster();
   out <<  ")";
}


tarch::parallel::messages::ActivationMessage::PersistentRecords tarch::parallel::messages::ActivationMessage::getPersistentRecords() const {
   return _persistentRecords;
}

tarch::parallel::messages::ActivationMessagePacked tarch::parallel::messages::ActivationMessage::convert() const{
   return ActivationMessagePacked(
      getNewMaster()
   );
}

#ifdef Parallel
   tarch::logging::Log tarch::parallel::messages::ActivationMessage::_log( "tarch::parallel::messages::ActivationMessage" );
   
   MPI_Datatype tarch::parallel::messages::ActivationMessage::Datatype = 0;
   MPI_Datatype tarch::parallel::messages::ActivationMessage::FullDatatype = 0;
   
   
   void tarch::parallel::messages::ActivationMessage::initDatatype() {
      {
         ActivationMessage dummyActivationMessage[2];
         
         const int Attributes = 2;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //newMaster
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //newMaster
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyActivationMessage[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyActivationMessage[0]._persistentRecords._newMaster))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyActivationMessage[1]._persistentRecords._newMaster))), 		&disp[1] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &ActivationMessage::Datatype );
         MPI_Type_commit( &ActivationMessage::Datatype );
         
      }
      {
         ActivationMessage dummyActivationMessage[2];
         
         const int Attributes = 2;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //newMaster
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //newMaster
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyActivationMessage[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyActivationMessage[0]._persistentRecords._newMaster))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyActivationMessage[1]._persistentRecords._newMaster))), 		&disp[1] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &ActivationMessage::FullDatatype );
         MPI_Type_commit( &ActivationMessage::FullDatatype );
         
      }
      
   }
   
   
   void tarch::parallel::messages::ActivationMessage::shutdownDatatype() {
      MPI_Type_free( &ActivationMessage::Datatype );
      MPI_Type_free( &ActivationMessage::FullDatatype );
      
   }
   
   void tarch::parallel::messages::ActivationMessage::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
         msg << "was not able to send message tarch::parallel::messages::ActivationMessage "
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
            msg << "testing for finished send task for tarch::parallel::messages::ActivationMessage "
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
            "tarch::parallel::messages::ActivationMessage",
            "send(int)", destination,tag,1
            );
            triggeredTimeoutWarning = true;
         }
         if (
            tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
            (clock()>timeOutShutdown)
         ) {
            tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "tarch::parallel::messages::ActivationMessage",
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
   
   
   
   void tarch::parallel::messages::ActivationMessage::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
         msg << "failed to start to receive tarch::parallel::messages::ActivationMessage from node "
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
            msg << "testing for finished receive task for tarch::parallel::messages::ActivationMessage failed: "
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
            "tarch::parallel::messages::ActivationMessage",
            "receive(int)", source,tag,1
            );
            triggeredTimeoutWarning = true;
         }
         if (
            tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
            (clock()>timeOutShutdown)
         ) {
            tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "tarch::parallel::messages::ActivationMessage",
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
   
   
   
   bool tarch::parallel::messages::ActivationMessage::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
   
   int tarch::parallel::messages::ActivationMessage::getSenderRank() const {
      assertion( _senderRank!=-1 );
      return _senderRank;
      
   }
#endif


tarch::parallel::messages::ActivationMessagePacked::PersistentRecords::PersistentRecords() {
   
}


tarch::parallel::messages::ActivationMessagePacked::PersistentRecords::PersistentRecords(const int& newMaster):
_newMaster(newMaster) {
   
}


 int tarch::parallel::messages::ActivationMessagePacked::PersistentRecords::getNewMaster() const  {
   return _newMaster;
}



 void tarch::parallel::messages::ActivationMessagePacked::PersistentRecords::setNewMaster(const int& newMaster)  {
   _newMaster = newMaster;
}


tarch::parallel::messages::ActivationMessagePacked::ActivationMessagePacked() {
   
}


tarch::parallel::messages::ActivationMessagePacked::ActivationMessagePacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._newMaster) {
   
}


tarch::parallel::messages::ActivationMessagePacked::ActivationMessagePacked(const int& newMaster):
_persistentRecords(newMaster) {
   
}


tarch::parallel::messages::ActivationMessagePacked::~ActivationMessagePacked() { }


 int tarch::parallel::messages::ActivationMessagePacked::getNewMaster() const  {
   return _persistentRecords._newMaster;
}



 void tarch::parallel::messages::ActivationMessagePacked::setNewMaster(const int& newMaster)  {
   _persistentRecords._newMaster = newMaster;
}




std::string tarch::parallel::messages::ActivationMessagePacked::toString() const {
   std::ostringstream stringstr;
   toString(stringstr);
   return stringstr.str();
}

void tarch::parallel::messages::ActivationMessagePacked::toString (std::ostream& out) const {
   out << "("; 
   out << "newMaster:" << getNewMaster();
   out <<  ")";
}


tarch::parallel::messages::ActivationMessagePacked::PersistentRecords tarch::parallel::messages::ActivationMessagePacked::getPersistentRecords() const {
   return _persistentRecords;
}

tarch::parallel::messages::ActivationMessage tarch::parallel::messages::ActivationMessagePacked::convert() const{
   return ActivationMessage(
      getNewMaster()
   );
}

#ifdef Parallel
   tarch::logging::Log tarch::parallel::messages::ActivationMessagePacked::_log( "tarch::parallel::messages::ActivationMessagePacked" );
   
   MPI_Datatype tarch::parallel::messages::ActivationMessagePacked::Datatype = 0;
   MPI_Datatype tarch::parallel::messages::ActivationMessagePacked::FullDatatype = 0;
   
   
   void tarch::parallel::messages::ActivationMessagePacked::initDatatype() {
      {
         ActivationMessagePacked dummyActivationMessagePacked[2];
         
         const int Attributes = 2;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //newMaster
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //newMaster
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyActivationMessagePacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyActivationMessagePacked[0]._persistentRecords._newMaster))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyActivationMessagePacked[1]._persistentRecords._newMaster))), 		&disp[1] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &ActivationMessagePacked::Datatype );
         MPI_Type_commit( &ActivationMessagePacked::Datatype );
         
      }
      {
         ActivationMessagePacked dummyActivationMessagePacked[2];
         
         const int Attributes = 2;
         MPI_Datatype subtypes[Attributes] = {
            MPI_INT,		 //newMaster
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            1,		 //newMaster
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyActivationMessagePacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyActivationMessagePacked[0]._persistentRecords._newMaster))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyActivationMessagePacked[1]._persistentRecords._newMaster))), 		&disp[1] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &ActivationMessagePacked::FullDatatype );
         MPI_Type_commit( &ActivationMessagePacked::FullDatatype );
         
      }
      
   }
   
   
   void tarch::parallel::messages::ActivationMessagePacked::shutdownDatatype() {
      MPI_Type_free( &ActivationMessagePacked::Datatype );
      MPI_Type_free( &ActivationMessagePacked::FullDatatype );
      
   }
   
   void tarch::parallel::messages::ActivationMessagePacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
         msg << "was not able to send message tarch::parallel::messages::ActivationMessagePacked "
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
            msg << "testing for finished send task for tarch::parallel::messages::ActivationMessagePacked "
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
            "tarch::parallel::messages::ActivationMessagePacked",
            "send(int)", destination,tag,1
            );
            triggeredTimeoutWarning = true;
         }
         if (
            tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
            (clock()>timeOutShutdown)
         ) {
            tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "tarch::parallel::messages::ActivationMessagePacked",
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
   
   
   
   void tarch::parallel::messages::ActivationMessagePacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
         msg << "failed to start to receive tarch::parallel::messages::ActivationMessagePacked from node "
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
            msg << "testing for finished receive task for tarch::parallel::messages::ActivationMessagePacked failed: "
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
            "tarch::parallel::messages::ActivationMessagePacked",
            "receive(int)", source,tag,1
            );
            triggeredTimeoutWarning = true;
         }
         if (
            tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
            (clock()>timeOutShutdown)
         ) {
            tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "tarch::parallel::messages::ActivationMessagePacked",
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
   
   
   
   bool tarch::parallel::messages::ActivationMessagePacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
   
   int tarch::parallel::messages::ActivationMessagePacked::getSenderRank() const {
      assertion( _senderRank!=-1 );
      return _senderRank;
      
   }
#endif



