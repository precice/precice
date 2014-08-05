#include "fsi/FSIDummyAImplementation.h"

#include <iostream>

fsi::FSIDummyAImplementation* fsi::FSIDummyAImplementation::singleton=NULL;
fsi::FSIDummyAImplementation::FSIDummyAImplementation():
		_localCoordIds(NULL),
		_numberOfPoints(0){
	pthread_mutex_init(&_mutex, NULL);
	std::cout<<"initializing FSIA"<<std::endl;
	singleton=this;
	//gatherMids();
	//gatherDomainDescriptions();
}

fsi::FSIDummyAImplementation::~FSIDummyAImplementation(){
	pthread_mutex_destroy(&_mutex);
}

//extern "C" void
//#ifdef _WIN32
//MAIN_LOOP(bool joinable);
//#else
//main_loop_(bool joinable);
//#endif

//int main(){
//#ifdef _WIN32
//	MAIN_LOOP(false);
//#else
//	main_loop_(false);
//#endif
//}

void fsi::FSIDummyAImplementation::setCoordIds(int* coordIds,const int numberOfPoints){
	pthread_mutex_lock(&_mutex);
	_localCoordIds=coordIds;
	_numberOfPoints=numberOfPoints;
	for(unsigned int i = 0 ; i< numberOfPoints ;i++)
		_global2LocalCoordMapping[coordIds[i]]=i;
	pthread_mutex_unlock(&_mutex);
}

void fsi::FSIDummyAImplementation::startDataTransfer(){

}
void fsi::FSIDummyAImplementation::endDataTransfer(int& ack){

}
void fsi::FSIDummyAImplementation::setData(double* data){
	_data=data;
}
void fsi::FSIDummyAImplementation::transferGlobalIds(){
	//double coordinates [2];
	//coordinates[0]=0.0;
	//coordinates[1]=0.0;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	if(rank==0){
		_b->transferCoordinatesParallel(
				&_globalCoordIds[0],
				_globalCoordIds.size(),
				&_globaOffsets[0],
				_globaOffsets.size(),
				&_mids[0],
				_mids.size()
		);
	}
}
void fsi::FSIDummyAImplementation::receiveAllData(){
	//pthread_mutex_lock(&_mutex);
	//double coordinates [2];
	//coordinates[0]=0.0;
	//coordinates[1]=0.0;
	int ack=-1;
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	if(rank==0){
		_b->startDataTransferParallel();
		_b->endDataTransfer(ack);
	}
	//pthread_mutex_unlock(&_mutex);
}

extern std::string retrieveSocketAddress();

void fsi::FSIDummyAImplementation::gatherDomainDescriptions(){
	int comm_size=0;
	int rank = 0;
	MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	if(comm_size>0)
		if(rank!=0){

			MPI_Send(&_numberOfPoints,1,MPI_INT,0,1000,MPI_COMM_WORLD);
			MPI_Send(&_localCoordIds,_numberOfPoints,MPI_INT,0,1001,MPI_COMM_WORLD);
		}else{
			//			int pieces=
			//					_parameters.parallel.numProcessors[0]*
			//					_parameters.parallel.numProcessors[1]*
			//					_parameters.parallel.numProcessors[2];
			int prefixSum=0;
			for(unsigned int i=0;i<_numberOfPoints;i++)
				_globalCoordIds.push_back(_localCoordIds[i]);
			_globaOffsets.push_back(prefixSum);
			//_globaOffsets.push_back(prefixSum);
			prefixSum+=_numberOfPoints;
			for(int i=1;i<comm_size;i++)
			{
				_globaOffsets.push_back(prefixSum);
				int numberOfPoints=0;
				MPI_Status status;
				MPI_Recv (&numberOfPoints,1, MPI_INT,i, 1000, MPI_COMM_WORLD,&status);
				//assert(number_of_bytes>0);
				std::vector<int> coordIds(_numberOfPoints);
				MPI_Recv (&coordIds[0],numberOfPoints,MPI_INT,i,1001,MPI_COMM_WORLD,&status);
				for(unsigned int i=0;i<numberOfPoints;i++)
					_globalCoordIds.push_back(coordIds[i]);
				prefixSum+=numberOfPoints;

			}

		}


	//	std::vector<int> startBuffer;
	//	std::vector<int> endBuffer;
	//	startBuffer.resize(3);
	//	endBuffer.resize(3);
	//	if(_parameters.parallel.rank==0){
	//		_startDomain.resize(3*_parameters.parallel.numProcessors[0]*_parameters.parallel.numProcessors[1]*_parameters.parallel.numProcessors[2]);
	//		_endDomain.resize(3*_parameters.parallel.numProcessors[0]*_parameters.parallel.numProcessors[1]*_parameters.parallel.numProcessors[2]);
	//	}
	//	startBuffer[0]=_parameters.parallel.firstCorner[0];
	//
	//	startBuffer[1]=_parameters.parallel.firstCorner[1];
	//
	//	startBuffer[2]=_parameters.parallel.firstCorner[2];
	//
	//	endBuffer[0]=_parameters.parallel.firstCorner[0]+_parameters.parallel.localSize[0]+1;
	//
	//	endBuffer[1]=_parameters.parallel.firstCorner[1]+_parameters.parallel.localSize[1]+1;
	//
	//	endBuffer[2]=_parameters.parallel.firstCorner[2]+_parameters.parallel.localSize[2]+1;
	//	MPI_Gather(&startBuffer[0],3,MPI_INT,&_startDomain[0],3,MPI_INT,0,MPI_COMM_WORLD);
	//	MPI_Gather(&endBuffer[0],3,MPI_INT,&_endDomain[0],3,MPI_INT,0,MPI_COMM_WORLD);

}
void fsi::FSIDummyAImplementation::gatherMids(){
	std::string mid=retrieveSocketAddress();
	_mids.push_back(mid);
	int comm_size=0;
	int rank = 0;
	MPI_Comm_size(MPI_COMM_WORLD,&comm_size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	if(comm_size>0)
		if(rank!=0){

			std::vector<char> bytes(mid.begin(), mid.end());
			bytes.push_back('\0');
			int num_of_bytes=(int)bytes.size();
			//assert(num_of_bytes>0);
			MPI_Send(&num_of_bytes,1,MPI_INT,0,1000,MPI_COMM_WORLD);
			MPI_Send(&bytes[0],num_of_bytes,MPI_BYTE,0,1000,MPI_COMM_WORLD);
		}else{
			//			int pieces=
			//					_parameters.parallel.numProcessors[0]*
			//					_parameters.parallel.numProcessors[1]*
			//					_parameters.parallel.numProcessors[2];
			for(int i=1;i<comm_size;i++)
			{
				int number_of_bytes=0;
				MPI_Status status;
				MPI_Recv (&number_of_bytes,1, MPI_INT,i, 1000, MPI_COMM_WORLD,&status);
				//assert(number_of_bytes>0);
				std::vector<char> buff(number_of_bytes);
				MPI_Recv (&buff[0],number_of_bytes,MPI_BYTE,i,1000,MPI_COMM_WORLD,&status);
				std::cout<<"rank:"<<i<<" mid:"<<std::string(&buff[0])<<std::endl;
				_mids.push_back(std::string(&buff[0]));
			}
		}

}
void fsi::FSIDummyAImplementation::dataAck(int& ack){
	pthread_mutex_lock(&_mutex);
	ack=1;
	pthread_mutex_unlock(&_mutex);
}

void fsi::FSIDummyAImplementation::transferData(
		const int* coordIds, const int coordIds_len,
		const double* data, const int data_len){
	pthread_mutex_lock(&_mutex);
	for(unsigned int i=0;i<data_len;i++)
		_data[_global2LocalCoordMapping[coordIds[i]]]=data[i];
	pthread_mutex_unlock(&_mutex);
}


