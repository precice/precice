#include <mpi.h>
#include "fsi/FSIDummyBImplementation.h"
#include "fsi/FSIDummyCommunicator.h"

#include <assert.h>
fsi::FSIDummyBImplementation* fsi::FSIDummyBImplementation::singleton=NULL;
fsi::FSIDummyBImplementation::FSIDummyBImplementation():
	_pointsSize(0),
	_initialized(false),
	_coordIds(NULL),
	_data(NULL){
	singleton=this;
	//	_localCoordinatesX(NULL),
//	_localCoordinatesY(NULL){
	pthread_mutex_init(&_mutex, NULL);
}

fsi::FSIDummyBImplementation::~FSIDummyBImplementation(){
	pthread_mutex_destroy(&_mutex);
}

//extern "C" void
//#ifdef _WIN32
//MAIN_LOOP(bool joinable);
//#else
//main_loop_(bool joinable);
//#endif
//
//int main(){
//#ifdef _WIN32
//	MAIN_LOOP(false);
//#else
//	main_loop_(false);
//#endif
//}

void fsi::FSIDummyBImplementation::setData(double* data){
	_data=data;
	_initialized=true;
}

void fsi::FSIDummyBImplementation::setCoordinates(
		int* coordIds,
//		double* coordinatesX,
//		double* coordinatesY,
		const int pointsSize){
	pthread_mutex_lock(&_mutex);
//	assert(coordinatesX);
//	assert(coordinatesY);
//	assert(_localIds);
	_coordIds=coordIds;
//	_localCoordinatesX=coordinatesX;
//	_localCoordinatesY=coordinatesY;
	_pointsSize=pointsSize;
	for(unsigned int i = 0 ; i< _pointsSize ;i++)
		_global2LocalCoordMapping[_coordIds[i]]=i;
	pthread_mutex_unlock(&_mutex);

}

void fsi::FSIDummyBImplementation::transferCoordinates(
		  const int* coordIds,
		  const int coord_len,
		  const int* offsets,
		  const int offsets_len,
		  const std::string* commids,
		  const int commids_len
){
	pthread_mutex_lock(&_mutex);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	int commId=-1;
	for(int i=0;i<coord_len;i++){
		for(int k=0;k<_pointsSize;k++)
			if(coordIds[i]==_coordIds[k]){
				commId=getHostId(i,offsets,offsets_len);

				assert(commId>=0);
				getCommunicator(commId,commids[commId])->setPointMapping(&_global2LocalCoordMapping);
				getCommunicator(commId,commids[commId])->addPointId(coordIds[i]);
				getCommunicator(commId,commids[commId])->setData(_data);
			}
		//std::cout<<"coord["<<i<<"]:"<<coord[i]<<std::endl;

	}


	//std::cout<<"end receiving coord"<<std::endl;
	pthread_mutex_unlock(&_mutex);
}

void fsi::FSIDummyBImplementation::flush(){
	pthread_mutex_lock(&_mutex);
	__gnu_cxx::hash_map<int,fsi::FSIDummyCommunicator*>::iterator it;
	for(it=_comms.begin();it!=_comms.end();++it){
		(*it).second->flush();
	}
	pthread_mutex_unlock(&_mutex);
}
fsi::FSIDummyCommunicator* fsi::FSIDummyBImplementation::getCommunicator(const int commId,const std::string commid){
	__gnu_cxx::hash_map<int,fsi::FSIDummyCommunicator*>::iterator it = _comms.find(commId);
	if(it==_comms.end()){
		_comms[commId]=new fsi::FSIDummyCommunicator(commid);
	}
	return _comms[commId];
}
const int fsi::FSIDummyBImplementation::getHostId(
		const int pointId,
		const int* offsets,
		const int offset_len) const {
	int res=0;
	//std::cout<<"offsets:"<<offsets[0]<<"point id:"<<pointId<<std::endl;
	for(unsigned int i=0;i<offset_len;i++){
		if(offsets[i]>pointId)
			return res-1;
		res++;
	}
	return res-1;
}

void fsi::FSIDummyBImplementation::startDataTransfer(){
	flush();
	MPI_Barrier(MPI_COMM_WORLD);
}
void fsi::FSIDummyBImplementation::endDataTransfer(int& ack){
	//fsi::FSIDummyBImplementation::MPI_Barrier();
	pthread_mutex_lock(&_mutex);
	ack=1;
	pthread_mutex_unlock(&_mutex);
}
void fsi::FSIDummyBImplementation::transferData(
		const int* coordIds,const int coordIds_len,
		const double* data,const int data_len){
	// @todo Insert your code here
}
