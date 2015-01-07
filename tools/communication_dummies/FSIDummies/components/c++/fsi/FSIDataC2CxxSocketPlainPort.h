#ifndef FSI_FSIDATAC2SOCKETPLAINPORT_H_
#define FSI_FSIDATAC2SOCKETPLAINPORT_H_ 

#include <map>
#include <string>


extern "C" {
#ifdef _WIN32
void FSI_FSIDATAC2SOCKET_PLAIN_PORT_CREATE_CLIENT_INSTANCE(long long* ptr,char* host,int& port,int& buffer_size);
void FSI_FSIDATAC2SOCKET_PLAIN_PORT_CREATE_SERVER_INSTANCE(long long* ptr,int& port,int& buffer_size);
void FSI_FSIDATAC2SOCKET_PLAIN_PORT_DESTROY_INSTANCE(long long* ptr);
#else
void fsi_fsidatac2socket_plain_port_create_client_instance_(long long* ptr,char* host,int& port,int& buffer_size);
void fsi_fsidatac2socket_plain_port_create_server_instance_(long long* ptr,int& port,int& buffer_size);
void fsi_fsidatac2socket_plain_port_destroy_instance_(long long* ptr);
#endif
#ifdef _WIN32
void FSI_FSIDATAC2SOCKET_PLAIN_PORT_TRANSFERDATA(long long* ref,int* coordId,int* coordId_len,double* data, int* data_len);

#else
void fsi_fsidatac2socket_plain_port_transferdata_(long long* ref,int* coordId,int* coordId_len,double* data, int* data_len);

#endif
#ifdef _WIN32
void FSI_FSIDATAC2SOCKET_PLAIN_PORT_DATAACK(long long* ref,int* ack);

#else
void fsi_fsidatac2socket_plain_port_dataack_(long long* ref,int* ack);

#endif
}
#endif
