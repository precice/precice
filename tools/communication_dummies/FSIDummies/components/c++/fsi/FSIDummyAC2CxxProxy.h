#ifndef __INCLUDE_GUARD_FULL_QUALIFIED_NAME__C2CXXPROXY_H_
#define __INCLUDE_GUARD_FULL_QUALIFIED_NAME__C2CXXPROXY_H_ 

#ifdef _WIN32
  #include <winsock2.h>
#endif

#ifdef Parallel
#include <mpi.h>
#endif
#include <pthread.h>
#include <string>
extern "C"{
#ifdef _WIN32
void SOCKET_CLIENT_LOOP();
#else
void socket_client_loop_();
#endif
 

#ifdef _WIN32
void SOCKET_SERVER_LOOP();
#else
void socket_server_loop_();
#endif

struct FSI_FSIDUMMYA_arg{
     void *ref;
     #ifdef _WIN32
     SOCKET
     #else
     int
     #endif
     java_serverfd;
     #ifdef _WIN32
     SOCKET
     #else
     int
     #endif
     java_clientfd;
     #ifdef _WIN32
     SOCKET
     #else
     int
     #endif
     daemon_serverfd;
     #ifdef _WIN32
     SOCKET
     #else
     int
     #endif
     daemon_clientfd;
     int buffer_size;
     std::string hostname;
     std::string client_port;
     std::string daemon_port;
     bool java_client_flag;
     bool joinable;
     int number_of_workers;
     const char* xml;
#ifdef Parallel
	 MPI_Comm communicator;
#endif     
};
}
#endif
