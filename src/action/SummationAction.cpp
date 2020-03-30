#include "SummationAction.hpp"
#include "mesh/Data.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Mesh.hpp"

namespace precice {
namespace action{

SummationAction::SummationAction(
	Timing				timing,
	std::vector<int>	sourceDataID,
	int					targetDataID,
	const mesh::PtrMesh	&mesh)
	: Action(timing, mesh, mapping::Mapping::MeshRequirement::VERTEX), _targetData(mesh->data(targetDataID))
{

	for(auto sourceID : sourceDataID){
		_sourceData.push_back(mesh->data(sourceID));
	}

	for(const auto & source : _sourceData) {
		PRECICE_ASSERT(source->getDimensions() == _targetData->getDimensions(), source->getDimensions(), _targetData->getDimensions());
	}

}

void SummationAction::performAction(
	double time,
	double dt,
	double computedPartFullDt,
	double fullDt)
{
	PRECICE_TRACE();
	auto & targetValues = _targetData->values();
	for(const auto & source : _sourceData){
		targetValues += source->values();
	}
}

}
}