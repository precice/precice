#include "fsi/FSIDummyCommunicator.h"

fsi::FSIDummyCommunicator::FSIDummyCommunicator(std::string hostname):
_hostname(hostname),
_data(NULL)
{

}


void fsi::FSIDummyCommunicator::setData(double* data){
	_data=data;
}
void fsi::FSIDummyCommunicator::addPointId(const int pointId){
	_pointIds.push_back(pointId);
}
void fsi::FSIDummyCommunicator::flush(){
	std::vector<double> data2Transfer;
	data2Transfer.resize(_pointIds.size());
	for(unsigned int i=0;i<_pointIds.size();i++)
	{
		data2Transfer[i]=_data[_pointIds[i]];
	}
	// do transfer

}
