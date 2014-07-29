#include "fsi/FSIDummyBImplementation.h"
#include "fsi/FSIDummyCommunicator.h"
#include <assert.h>
fsi::FSIDummyBImplementation::FSIDummyBImplementation():
	_pointsSize(0),
	_initialized(false),
	_localIds(NULL),
	_data(NULL){
	//	_localCoordinatesX(NULL),
//	_localCoordinatesY(NULL){

}

fsi::FSIDummyBImplementation::~FSIDummyBImplementation(){

}

extern "C" void
#ifdef _WIN32
MAIN_LOOP(bool joinable);
#else
main_loop_(bool joinable);
#endif

int main(){
#ifdef _WIN32
	MAIN_LOOP(false);
#else
	main_loop_(false);
#endif
}

void fsi::FSIDummyBImplementation::setData(double* data){
	_data=data;
	_initialized=true;
}

void fsi::FSIDummyBImplementation::setCoordinates(
		int* localIds,
//		double* coordinatesX,
//		double* coordinatesY,
		const int pointsSize){
//	assert(coordinatesX);
//	assert(coordinatesY);
//	assert(_localIds);
	_localIds=localIds;
//	_localCoordinatesX=coordinatesX;
//	_localCoordinatesY=coordinatesY;
	_pointsSize=pointsSize;

}

void fsi::FSIDummyBImplementation::transferCoordinates(
		  const int* coordIds,
		  const int coord_len,
		  const int* offsets,
		  const int offsets_len,
		  const std::string* commids,
		  const int commids_len
){
	std::cout<<"start receiving coord"<<std::endl;
	int commId=-1;
	for(int i=0;i<coord_len;i++){
		for(int k=0;k<_pointsSize;k++)
			if(coordIds[i]==_localIds[k]){
				commId=getHostId(i,offsets,offsets_len);

				assert(commId>=0);

				getCommunicator(commId,commids[commId])->addPointId(coordIds[i]);
			}
		//std::cout<<"coord["<<i<<"]:"<<coord[i]<<std::endl;

	}
	//std::cout<<"end receiving coord"<<std::endl;

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
	for(unsigned int i=0;i<offset_len;i++){
		if(offsets[i]>pointId)
			return res;
		res++;
	}
	return -1;
}
void fsi::FSIDummyBImplementation::transferData(const double* data, const int data_len){
	// @todo Insert your code here
}
