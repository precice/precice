#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/parallel/messages/RegisterAtNodePoolMessage.h"

tarch::parallel::messages::RegisterAtNodePoolMessage::PersistentRecords::PersistentRecords() {
   
}


tarch::parallel::messages::RegisterAtNodePoolMessage::PersistentRecords::PersistentRecords(const tarch::la::Vector<MPI_MAX_NAME_STRING_ADDED_ONE,short int>& nodeName):
_nodeName(nodeName) {
   
}


 tarch::la::Vector<MPI_MAX_NAME_STRING_ADDED_ONE,short int> tarch::parallel::messages::RegisterAtNodePoolMessage::PersistentRecords::getNodeName() const  {
   return _nodeName;
}



 void tarch::parallel::messages::RegisterAtNodePoolMessage::PersistentRecords::setNodeName(const tarch::la::Vector<MPI_MAX_NAME_STRING_ADDED_ONE,short int>& nodeName)  {
   _nodeName = (nodeName);
}


tarch::parallel::messages::RegisterAtNodePoolMessage::RegisterAtNodePoolMessage() {
   
}


tarch::parallel::messages::RegisterAtNodePoolMessage::RegisterAtNodePoolMessage(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._nodeName) {
   
}


tarch::parallel::messages::RegisterAtNodePoolMessage::RegisterAtNodePoolMessage(const tarch::la::Vector<MPI_MAX_NAME_STRING_ADDED_ONE,short int>& nodeName):
_persistentRecords(nodeName) {
   
}


tarch::parallel::messages::RegisterAtNodePoolMessage::~RegisterAtNodePoolMessage() { }


 tarch::la::Vector<MPI_MAX_NAME_STRING_ADDED_ONE,short int> tarch::parallel::messages::RegisterAtNodePoolMessage::getNodeName() const  {
   return _persistentRecords._nodeName;
}



 void tarch::parallel::messages::RegisterAtNodePoolMessage::setNodeName(const tarch::la::Vector<MPI_MAX_NAME_STRING_ADDED_ONE,short int>& nodeName)  {
   _persistentRecords._nodeName = (nodeName);
}



 short int tarch::parallel::messages::RegisterAtNodePoolMessage::getNodeName(int elementIndex) const  {
   assertion(elementIndex>=0);
   assertion(elementIndex<MPI_MAX_NAME_STRING_ADDED_ONE);
   return _persistentRecords._nodeName[elementIndex];
   
}



 void tarch::parallel::messages::RegisterAtNodePoolMessage::setNodeName(int elementIndex, const short int& nodeName)  {
   assertion(elementIndex>=0);
   assertion(elementIndex<MPI_MAX_NAME_STRING_ADDED_ONE);
   _persistentRecords._nodeName[elementIndex]= nodeName;
   
}




std::string tarch::parallel::messages::RegisterAtNodePoolMessage::toString() const {
   std::ostringstream stringstr;
   toString(stringstr);
   return stringstr.str();
}

void tarch::parallel::messages::RegisterAtNodePoolMessage::toString (std::ostream& out) const {
   out << "("; 
   out << "nodeName:[";
   for (int i = 0; i < MPI_MAX_NAME_STRING_ADDED_ONE-1; i++) {
      out << getNodeName(i) << ",";
   }
   out << getNodeName(MPI_MAX_NAME_STRING_ADDED_ONE-1) << "]";
   out <<  ")";
}


tarch::parallel::messages::RegisterAtNodePoolMessage::PersistentRecords tarch::parallel::messages::RegisterAtNodePoolMessage::getPersistentRecords() const {
   return _persistentRecords;
}

tarch::parallel::messages::RegisterAtNodePoolMessagePacked tarch::parallel::messages::RegisterAtNodePoolMessage::convert() const{
   return RegisterAtNodePoolMessagePacked(
      getNodeName()
   );
}

#ifdef Parallel
   tarch::logging::Log tarch::parallel::messages::RegisterAtNodePoolMessage::_log( "tarch::parallel::messages::RegisterAtNodePoolMessage" );
   
   MPI_Datatype tarch::parallel::messages::RegisterAtNodePoolMessage::Datatype = 0;
   MPI_Datatype tarch::parallel::messages::RegisterAtNodePoolMessage::FullDatatype = 0;
   
   
   void tarch::parallel::messages::RegisterAtNodePoolMessage::initDatatype() {
      {
         RegisterAtNodePoolMessage dummyRegisterAtNodePoolMessage[2];
         
         const int Attributes = 2;
         MPI_Datatype subtypes[Attributes] = {
            MPI_SHORT,		 //nodeName
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            MPI_MAX_NAME_STRING_ADDED_ONE,		 //nodeName
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRegisterAtNodePoolMessage[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRegisterAtNodePoolMessage[0]._persistentRecords._nodeName[0]))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&dummyRegisterAtNodePoolMessage[1]._persistentRecords._nodeName[0])), 		&disp[1] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &RegisterAtNodePoolMessage::Datatype );
         MPI_Type_commit( &RegisterAtNodePoolMessage::Datatype );
         
      }
      {
         RegisterAtNodePoolMessage dummyRegisterAtNodePoolMessage[2];
         
         const int Attributes = 2;
         MPI_Datatype subtypes[Attributes] = {
            MPI_SHORT,		 //nodeName
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            MPI_MAX_NAME_STRING_ADDED_ONE,		 //nodeName
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRegisterAtNodePoolMessage[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRegisterAtNodePoolMessage[0]._persistentRecords._nodeName[0]))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&dummyRegisterAtNodePoolMessage[1]._persistentRecords._nodeName[0])), 		&disp[1] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &RegisterAtNodePoolMessage::FullDatatype );
         MPI_Type_commit( &RegisterAtNodePoolMessage::FullDatatype );
         
      }
      
   }
   
   
   void tarch::parallel::messages::RegisterAtNodePoolMessage::shutdownDatatype() {
      MPI_Type_free( &RegisterAtNodePoolMessage::Datatype );
      MPI_Type_free( &RegisterAtNodePoolMessage::FullDatatype );
      
   }
   
   void tarch::parallel::messages::RegisterAtNodePoolMessage::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
         msg << "was not able to send message tarch::parallel::messages::RegisterAtNodePoolMessage "
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
            msg << "testing for finished send task for tarch::parallel::messages::RegisterAtNodePoolMessage "
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
            "tarch::parallel::messages::RegisterAtNodePoolMessage",
            "send(int)", destination,tag,1
            );
            triggeredTimeoutWarning = true;
         }
         if (
            tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
            (clock()>timeOutShutdown)
         ) {
            tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "tarch::parallel::messages::RegisterAtNodePoolMessage",
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
   
   
   
   void tarch::parallel::messages::RegisterAtNodePoolMessage::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
         msg << "failed to start to receive tarch::parallel::messages::RegisterAtNodePoolMessage from node "
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
            msg << "testing for finished receive task for tarch::parallel::messages::RegisterAtNodePoolMessage failed: "
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
            "tarch::parallel::messages::RegisterAtNodePoolMessage",
            "receive(int)", source,tag,1
            );
            triggeredTimeoutWarning = true;
         }
         if (
            tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
            (clock()>timeOutShutdown)
         ) {
            tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "tarch::parallel::messages::RegisterAtNodePoolMessage",
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
   
   
   
   bool tarch::parallel::messages::RegisterAtNodePoolMessage::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
   
   int tarch::parallel::messages::RegisterAtNodePoolMessage::getSenderRank() const {
      assertion( _senderRank!=-1 );
      return _senderRank;
      
   }
#endif


tarch::parallel::messages::RegisterAtNodePoolMessagePacked::PersistentRecords::PersistentRecords() {
   
}


tarch::parallel::messages::RegisterAtNodePoolMessagePacked::PersistentRecords::PersistentRecords(const tarch::la::Vector<MPI_MAX_NAME_STRING_ADDED_ONE,short int>& nodeName):
_nodeName(nodeName) {
   
}


 tarch::la::Vector<MPI_MAX_NAME_STRING_ADDED_ONE,short int> tarch::parallel::messages::RegisterAtNodePoolMessagePacked::PersistentRecords::getNodeName() const  {
   return _nodeName;
}



 void tarch::parallel::messages::RegisterAtNodePoolMessagePacked::PersistentRecords::setNodeName(const tarch::la::Vector<MPI_MAX_NAME_STRING_ADDED_ONE,short int>& nodeName)  {
   _nodeName = (nodeName);
}


tarch::parallel::messages::RegisterAtNodePoolMessagePacked::RegisterAtNodePoolMessagePacked() {
   
}


tarch::parallel::messages::RegisterAtNodePoolMessagePacked::RegisterAtNodePoolMessagePacked(const PersistentRecords& persistentRecords):
_persistentRecords(persistentRecords._nodeName) {
   
}


tarch::parallel::messages::RegisterAtNodePoolMessagePacked::RegisterAtNodePoolMessagePacked(const tarch::la::Vector<MPI_MAX_NAME_STRING_ADDED_ONE,short int>& nodeName):
_persistentRecords(nodeName) {
   
}


tarch::parallel::messages::RegisterAtNodePoolMessagePacked::~RegisterAtNodePoolMessagePacked() { }


 tarch::la::Vector<MPI_MAX_NAME_STRING_ADDED_ONE,short int> tarch::parallel::messages::RegisterAtNodePoolMessagePacked::getNodeName() const  {
   return _persistentRecords._nodeName;
}



 void tarch::parallel::messages::RegisterAtNodePoolMessagePacked::setNodeName(const tarch::la::Vector<MPI_MAX_NAME_STRING_ADDED_ONE,short int>& nodeName)  {
   _persistentRecords._nodeName = (nodeName);
}



 short int tarch::parallel::messages::RegisterAtNodePoolMessagePacked::getNodeName(int elementIndex) const  {
   assertion(elementIndex>=0);
   assertion(elementIndex<MPI_MAX_NAME_STRING_ADDED_ONE);
   return _persistentRecords._nodeName[elementIndex];
   
}



 void tarch::parallel::messages::RegisterAtNodePoolMessagePacked::setNodeName(int elementIndex, const short int& nodeName)  {
   assertion(elementIndex>=0);
   assertion(elementIndex<MPI_MAX_NAME_STRING_ADDED_ONE);
   _persistentRecords._nodeName[elementIndex]= nodeName;
   
}




std::string tarch::parallel::messages::RegisterAtNodePoolMessagePacked::toString() const {
   std::ostringstream stringstr;
   toString(stringstr);
   return stringstr.str();
}

void tarch::parallel::messages::RegisterAtNodePoolMessagePacked::toString (std::ostream& out) const {
   out << "("; 
   out << "nodeName:[";
   for (int i = 0; i < MPI_MAX_NAME_STRING_ADDED_ONE-1; i++) {
      out << getNodeName(i) << ",";
   }
   out << getNodeName(MPI_MAX_NAME_STRING_ADDED_ONE-1) << "]";
   out <<  ")";
}


tarch::parallel::messages::RegisterAtNodePoolMessagePacked::PersistentRecords tarch::parallel::messages::RegisterAtNodePoolMessagePacked::getPersistentRecords() const {
   return _persistentRecords;
}

tarch::parallel::messages::RegisterAtNodePoolMessage tarch::parallel::messages::RegisterAtNodePoolMessagePacked::convert() const{
   return RegisterAtNodePoolMessage(
      getNodeName()
   );
}

#ifdef Parallel
   tarch::logging::Log tarch::parallel::messages::RegisterAtNodePoolMessagePacked::_log( "tarch::parallel::messages::RegisterAtNodePoolMessagePacked" );
   
   MPI_Datatype tarch::parallel::messages::RegisterAtNodePoolMessagePacked::Datatype = 0;
   MPI_Datatype tarch::parallel::messages::RegisterAtNodePoolMessagePacked::FullDatatype = 0;
   
   
   void tarch::parallel::messages::RegisterAtNodePoolMessagePacked::initDatatype() {
      {
         RegisterAtNodePoolMessagePacked dummyRegisterAtNodePoolMessagePacked[2];
         
         const int Attributes = 2;
         MPI_Datatype subtypes[Attributes] = {
            MPI_SHORT,		 //nodeName
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            MPI_MAX_NAME_STRING_ADDED_ONE,		 //nodeName
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRegisterAtNodePoolMessagePacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRegisterAtNodePoolMessagePacked[0]._persistentRecords._nodeName[0]))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&dummyRegisterAtNodePoolMessagePacked[1]._persistentRecords._nodeName[0])), 		&disp[1] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &RegisterAtNodePoolMessagePacked::Datatype );
         MPI_Type_commit( &RegisterAtNodePoolMessagePacked::Datatype );
         
      }
      {
         RegisterAtNodePoolMessagePacked dummyRegisterAtNodePoolMessagePacked[2];
         
         const int Attributes = 2;
         MPI_Datatype subtypes[Attributes] = {
            MPI_SHORT,		 //nodeName
            MPI_UB		 // end/displacement flag
         };
         
         int blocklen[Attributes] = {
            MPI_MAX_NAME_STRING_ADDED_ONE,		 //nodeName
            1		 // end/displacement flag
         };
         
         MPI_Aint     disp[Attributes];
         
         MPI_Aint base;
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRegisterAtNodePoolMessagePacked[0]))), &base);
         MPI_Address( const_cast<void*>(static_cast<const void*>(&(dummyRegisterAtNodePoolMessagePacked[0]._persistentRecords._nodeName[0]))), 		&disp[0] );
         MPI_Address( const_cast<void*>(static_cast<const void*>(&dummyRegisterAtNodePoolMessagePacked[1]._persistentRecords._nodeName[0])), 		&disp[1] );
         
         for (int i=1; i<Attributes; i++) {
            assertion1( disp[i] > disp[i-1], i );
         }
         for (int i=0; i<Attributes; i++) {
            disp[i] -= base;
         }
         MPI_Type_struct( Attributes, blocklen, disp, subtypes, &RegisterAtNodePoolMessagePacked::FullDatatype );
         MPI_Type_commit( &RegisterAtNodePoolMessagePacked::FullDatatype );
         
      }
      
   }
   
   
   void tarch::parallel::messages::RegisterAtNodePoolMessagePacked::shutdownDatatype() {
      MPI_Type_free( &RegisterAtNodePoolMessagePacked::Datatype );
      MPI_Type_free( &RegisterAtNodePoolMessagePacked::FullDatatype );
      
   }
   
   void tarch::parallel::messages::RegisterAtNodePoolMessagePacked::send(int destination, int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
         msg << "was not able to send message tarch::parallel::messages::RegisterAtNodePoolMessagePacked "
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
            msg << "testing for finished send task for tarch::parallel::messages::RegisterAtNodePoolMessagePacked "
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
            "tarch::parallel::messages::RegisterAtNodePoolMessagePacked",
            "send(int)", destination,tag,1
            );
            triggeredTimeoutWarning = true;
         }
         if (
            tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
            (clock()>timeOutShutdown)
         ) {
            tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "tarch::parallel::messages::RegisterAtNodePoolMessagePacked",
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
   
   
   
   void tarch::parallel::messages::RegisterAtNodePoolMessagePacked::receive(int source, int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
         msg << "failed to start to receive tarch::parallel::messages::RegisterAtNodePoolMessagePacked from node "
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
            msg << "testing for finished receive task for tarch::parallel::messages::RegisterAtNodePoolMessagePacked failed: "
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
            "tarch::parallel::messages::RegisterAtNodePoolMessagePacked",
            "receive(int)", source,tag,1
            );
            triggeredTimeoutWarning = true;
         }
         if (
            tarch::parallel::Node::getInstance().isTimeOutDeadlockEnabled() &&
            (clock()>timeOutShutdown)
         ) {
            tarch::parallel::Node::getInstance().triggerDeadlockTimeOut(
            "tarch::parallel::messages::RegisterAtNodePoolMessagePacked",
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
   
   
   
   bool tarch::parallel::messages::RegisterAtNodePoolMessagePacked::isMessageInQueue(int tag, bool exchangeOnlyAttributesMarkedWithParallelise) {
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
   
   int tarch::parallel::messages::RegisterAtNodePoolMessagePacked::getSenderRank() const {
      assertion( _senderRank!=-1 );
      return _senderRank;
      
   }
#endif



