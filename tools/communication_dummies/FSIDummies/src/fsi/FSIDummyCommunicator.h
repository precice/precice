#ifndef FSI_DUMMY_COMMUNICATOR_H
#define FSI_DUMMY_COMMUNICATOR_H

#include <vector>
#include <string>
namespace fsi{
	class FSIDummyCommunicator;
}

class fsi::FSIDummyCommunicator{
private:
	std::vector<int> _pointIds;
	std::string _hostname;
	double* _data;
public:
	FSIDummyCommunicator(std::string hostname);
	void setData(double* data);
	void addPointId(const int pointId);
	void flush();
};

#endif
