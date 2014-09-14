#ifndef PRECICE_INITIALIZERC2SOCKETPLAINPORT_H_
#define PRECICE_INITIALIZERC2SOCKETPLAINPORT_H_ 

#include <map>
#include <string>


extern "C" {
#ifdef _WIN32
void PRECICE_INITIALIZERC2SOCKET_PLAIN_PORT_CREATE_CLIENT_INSTANCE(long long* ptr,char* host,int& port,int& buffer_size);
void PRECICE_INITIALIZERC2SOCKET_PLAIN_PORT_CREATE_SERVER_INSTANCE(long long* ptr,int& port,int& buffer_size);
void PRECICE_INITIALIZERC2SOCKET_PLAIN_PORT_DESTROY_INSTANCE(long long* ptr);
#else
void precice_initializerc2socket_plain_port_create_client_instance_(long long* ptr,char* host,int& port,int& buffer_size);
void precice_initializerc2socket_plain_port_create_server_instance_(long long* ptr,int& port,int& buffer_size);
void precice_initializerc2socket_plain_port_destroy_instance_(long long* ptr);
#endif
#ifdef _WIN32
void PRECICE_INITIALIZERC2SOCKET_PLAIN_PORT_INITIALIZEADDRESSES(long long* ref,char** addresses,int* addresses_len);

#else
void precice_initializerc2socket_plain_port_initializeaddresses_(long long* ref,char** addresses,int* addresses_len);

#endif
#ifdef _WIN32
void PRECICE_INITIALIZERC2SOCKET_PLAIN_PORT_INITIALIZEVERTEXES(long long* ref,int* vertexes,int* vertexes_len);

#else
void precice_initializerc2socket_plain_port_initializevertexes_(long long* ref,int* vertexes,int* vertexes_len);

#endif
}
#endif
