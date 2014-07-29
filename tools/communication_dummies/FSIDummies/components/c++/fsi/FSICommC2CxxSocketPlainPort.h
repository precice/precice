#ifndef FSI_FSICOMMC2SOCKETPLAINPORT_H_
#define FSI_FSICOMMC2SOCKETPLAINPORT_H_ 

#include <map>
#include <string>


extern "C" {
#ifdef _WIN32
void FSI_FSICOMMC2SOCKET_PLAIN_PORT_CREATE_CLIENT_INSTANCE(long long* ptr,char* host,int& port,int& buffer_size);
void FSI_FSICOMMC2SOCKET_PLAIN_PORT_CREATE_SERVER_INSTANCE(long long* ptr,int& port,int& buffer_size);
void FSI_FSICOMMC2SOCKET_PLAIN_PORT_DESTROY_INSTANCE(long long* ptr);
#else
void fsi_fsicommc2socket_plain_port_create_client_instance_(long long* ptr,char* host,int& port,int& buffer_size);
void fsi_fsicommc2socket_plain_port_create_server_instance_(long long* ptr,int& port,int& buffer_size);
void fsi_fsicommc2socket_plain_port_destroy_instance_(long long* ptr);
#endif
#ifdef _WIN32
void FSI_FSICOMMC2SOCKET_PLAIN_PORT_TRANSFERCOORDINATES(long long* ref,double* coord, int* coord_len);

#else
void fsi_fsicommc2socket_plain_port_transfercoordinates_(long long* ref,double* coord, int* coord_len);

#endif
#ifdef _WIN32
void FSI_FSICOMMC2SOCKET_PLAIN_PORT_TRANSFERDATA(long long* ref,double* data, int* data_len);

#else
void fsi_fsicommc2socket_plain_port_transferdata_(long long* ref,double* data, int* data_len);

#endif
}
#endif
