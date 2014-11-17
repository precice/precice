#ifndef PRECICE_MAINC2SOCKETPLAINPORT_H_
#define PRECICE_MAINC2SOCKETPLAINPORT_H_ 

#include <map>
#include <string>


extern "C" {
#ifdef _WIN32
void PRECICE_MAINC2SOCKET_PLAIN_PORT_CREATE_CLIENT_INSTANCE(long long* ptr,char* host,int& port,int& buffer_size);
void PRECICE_MAINC2SOCKET_PLAIN_PORT_CREATE_SERVER_INSTANCE(long long* ptr,int& port,int& buffer_size);
void PRECICE_MAINC2SOCKET_PLAIN_PORT_DESTROY_INSTANCE(long long* ptr);
#else
void precice_mainc2socket_plain_port_create_client_instance_(long long* ptr,char* host,int& port,int& buffer_size);
void precice_mainc2socket_plain_port_create_server_instance_(long long* ptr,int& port,int& buffer_size);
void precice_mainc2socket_plain_port_destroy_instance_(long long* ptr);
#endif
#ifdef _WIN32
void PRECICE_MAINC2SOCKET_PLAIN_PORT_MAIN(long long* ref);

#else
void precice_mainc2socket_plain_port_main_(long long* ref);

#endif
}
#endif
