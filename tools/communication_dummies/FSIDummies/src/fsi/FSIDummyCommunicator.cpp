#include "fsi/FSIDummyCommunicator.h"
#include "fsi/FSIDataCxx2SocketPlainPort.h"

#include <stdlib.h>
fsi::FSIDummyCommunicator::FSIDummyCommunicator(std::string hostname):
_hostname(hostname),
_data(NULL),
_comm(NULL),
_globalToLocalPointMapping(NULL)
{
}

fsi::FSIDummyCommunicator::~FSIDummyCommunicator(){
	if(_comm){
		disconnect();
	}
}

void fsi::FSIDummyCommunicator::connect(){
	if(_comm==NULL){
		std::string hostname = _hostname.substr(0,_hostname.find(":"));
		std::string port = _hostname.substr(_hostname.find(":")+1,_hostname.size()-_hostname.find(":"));
		const char* buffer_size = getenv("FSI_FSIDUMMYB_BUFFER_SIZE");
		_comm= new fsi::FSIDataCxx2SocketPlainPort(
				(char*)hostname.c_str(),atoi(port.c_str()),atoi(buffer_size)
		);
	}
}

void fsi::FSIDummyCommunicator::disconnect(){
	delete _comm;
	_comm = NULL;
}

void fsi::FSIDummyCommunicator::setData(double* data){
	_data=data;
}

void fsi::FSIDummyCommunicator::addPointId(const int pointId){
	_pointIds.push_back(pointId);
}

void fsi::FSIDummyCommunicator::setPointMapping(__gnu_cxx::hash_map<int,int>* globalToLocalPointMapping){
	_globalToLocalPointMapping= globalToLocalPointMapping;
}
void fsi::FSIDummyCommunicator::flush(){
	connect();
	_data2Transfer.resize(_pointIds.size());
	std::cout<<"data to send:"<<std::endl;
	for(unsigned int i=0;i<_pointIds.size();i++)
	{
		_data2Transfer[i]=_data[(*_globalToLocalPointMapping)[_pointIds[i]]];
		std::cout<<"data "<<i<<":"<<_data2Transfer[i]<<std::endl;
	}
	_comm->transferData(&_pointIds[0],_pointIds.size(),&_data2Transfer[0],_data2Transfer.size());
	int ack;
	_comm->dataAck(ack);
	std::cout<<"data flushed"<<std::endl;
	// do transfer

}
