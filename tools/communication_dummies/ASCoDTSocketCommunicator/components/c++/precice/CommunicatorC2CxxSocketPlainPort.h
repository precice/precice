#ifndef PRECICE_COMMUNICATORC2SOCKETPLAINPORT_H_
#define PRECICE_COMMUNICATORC2SOCKETPLAINPORT_H_ 

#include <map>
#include <string>


extern "C" {
#ifdef _WIN32
void PRECICE_COMMUNICATORC2SOCKET_PLAIN_PORT_CREATE_CLIENT_INSTANCE(long long* ptr,char* host,int& port,int& buffer_size);
void PRECICE_COMMUNICATORC2SOCKET_PLAIN_PORT_CREATE_SERVER_INSTANCE(long long* ptr,int& port,int& buffer_size);
void PRECICE_COMMUNICATORC2SOCKET_PLAIN_PORT_DESTROY_INSTANCE(long long* ptr);
#else
void precice_communicatorc2socket_plain_port_create_client_instance_(long long* ptr,char* host,int& port,int& buffer_size);
void precice_communicatorc2socket_plain_port_create_server_instance_(long long* ptr,int& port,int& buffer_size);
void precice_communicatorc2socket_plain_port_destroy_instance_(long long* ptr);
#endif
#ifdef _WIN32
void PRECICE_COMMUNICATORC2SOCKET_PLAIN_PORT_SETDATA(long long* ref,double* data,int* index,int* rank,int* tag);

#else
void precice_communicatorc2socket_plain_port_setdata_(long long* ref,double* data,int* index,int* rank,int* tag);

#endif
}
#endif
