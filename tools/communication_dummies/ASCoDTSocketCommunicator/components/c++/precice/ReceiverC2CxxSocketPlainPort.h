#ifndef PRECICE_RECEIVERC2SOCKETPLAINPORT_H_
#define PRECICE_RECEIVERC2SOCKETPLAINPORT_H_ 

#include <map>
#include <string>


extern "C" {
#ifdef _WIN32
void PRECICE_RECEIVERC2SOCKET_PLAIN_PORT_CREATE_CLIENT_INSTANCE(long long* ptr,char* host,int& port,int& buffer_size);
void PRECICE_RECEIVERC2SOCKET_PLAIN_PORT_CREATE_SERVER_INSTANCE(long long* ptr,int& port,int& buffer_size);
void PRECICE_RECEIVERC2SOCKET_PLAIN_PORT_DESTROY_INSTANCE(long long* ptr);
#else
void precice_receiverc2socket_plain_port_create_client_instance_(long long* ptr,char* host,int& port,int& buffer_size);
void precice_receiverc2socket_plain_port_create_server_instance_(long long* ptr,int& port,int& buffer_size);
void precice_receiverc2socket_plain_port_destroy_instance_(long long* ptr);
#endif
#ifdef _WIN32
void PRECICE_RECEIVERC2SOCKET_PLAIN_PORT_RECEIVE(long long* ref,double* data,int* index,int* rank,int* tag);

#else
void precice_receiverc2socket_plain_port_receive_(long long* ref,double* data,int* index,int* rank,int* tag);

#endif
}
#endif
