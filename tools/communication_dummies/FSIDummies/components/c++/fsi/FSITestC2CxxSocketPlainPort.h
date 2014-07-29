#ifndef FSI_FSITESTC2SOCKETPLAINPORT_H_
#define FSI_FSITESTC2SOCKETPLAINPORT_H_ 

#include <map>
#include <string>


extern "C" {
#ifdef _WIN32
void FSI_FSITESTC2SOCKET_PLAIN_PORT_CREATE_CLIENT_INSTANCE(long long* ptr,char* host,int& port,int& buffer_size);
void FSI_FSITESTC2SOCKET_PLAIN_PORT_CREATE_SERVER_INSTANCE(long long* ptr,int& port,int& buffer_size);
void FSI_FSITESTC2SOCKET_PLAIN_PORT_DESTROY_INSTANCE(long long* ptr);
#else
void fsi_fsitestc2socket_plain_port_create_client_instance_(long long* ptr,char* host,int& port,int& buffer_size);
void fsi_fsitestc2socket_plain_port_create_server_instance_(long long* ptr,int& port,int& buffer_size);
void fsi_fsitestc2socket_plain_port_destroy_instance_(long long* ptr);
#endif
#ifdef _WIN32
void FSI_FSITESTC2SOCKET_PLAIN_PORT_TEST(long long* ref);

#else
void fsi_fsitestc2socket_plain_port_test_(long long* ref);

#endif
}
#endif
