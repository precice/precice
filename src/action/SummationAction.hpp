#pragma once

#include "Action.hpp"
#include "logging/Logger.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {
namespace mesh {
class Edge;
class Triangle;
} // namespace mesh
} // namespace precice

namespace precice {
namespace action {

class SummationAction : public Action {
public:
	/**
	 * @brief Constructor
	 * 
	 * @param[in] timing When to apply the action
	 * @param[in] source1 First source data
	 * @param[in] source2 Second source data
	 * @param[in] target Data in which the action will be applied
	 * 
	 */
	SummationAction(
		Timing				timing,
		std::vector<int>	sourceDataID,
		int					targetDataID,
		const mesh::PtrMesh	&mesh);
	
	virtual ~SummationAction() {}
	 
	/**
	 * @brief Adds or substracts the source data to target data
	 * 
	 */
	virtual void performAction(
		double time,
		double dt,
		double computedPartFullDt,
		double fullDt);

private:
	logging::Logger _log{"action::SummationAction"};

	mesh::PtrData _targetData;
	std::vector<mesh::PtrData> _sourceData;

};

} // namespace precice
} // namespace action