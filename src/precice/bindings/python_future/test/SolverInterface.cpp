#include "precice/SolverInterface.hpp"

std::vector<double> fake_read_write_buffer;
int fake_dimensions;
int fake_mesh_id;

namespace precice {

namespace impl{
class SolverInterfaceImpl{};
}

SolverInterface:: SolverInterface
(
  const std::string& participantName,
  int                solverProcessIndex,
  int                solverProcessSize )
{
  fake_read_write_buffer = std::vector<double>();
  fake_dimensions = 3;
  fake_mesh_id = 0;
}

SolverInterface::~SolverInterface() = default;

void SolverInterface:: configure
(
  const std::string& configurationFileName )
{}

double SolverInterface:: initialize(){return -1;}

void SolverInterface:: initializeData()
{}

double SolverInterface:: advance
(
  double computedTimestepLength )
{return -1;}

void SolverInterface:: finalize()
{}

int SolverInterface:: getDimensions() const
{return fake_dimensions;}

bool SolverInterface:: isCouplingOngoing() const
{
  return 0;
}

bool SolverInterface:: isReadDataAvailable() const
{
  return 0;
}

bool SolverInterface:: isWriteDataRequired
(
  double computedTimestepLength ) const
{
  return 0;
}

bool SolverInterface:: isTimestepComplete() const
{
  return 0;
}

bool SolverInterface:: isActionRequired
(
  const std::string& action ) const
{
  return 0;
}

void SolverInterface:: fulfilledAction
(
  const std::string& action )
{}

bool SolverInterface:: hasMesh
(
  const std::string& meshName ) const
{
  return 0;
}

int SolverInterface:: getMeshID
(
  const std::string& meshName ) const
{
  return fake_mesh_id;
}

std::set<int> SolverInterface:: getMeshIDs() const
{
  return std::set<int>();
}

bool SolverInterface:: hasData
(
  const std::string& dataName, int meshID ) const
{
  return 0;
}

int SolverInterface:: getDataID
(
  const std::string& dataName, int meshID ) const
{
  return -1;
}

bool SolverInterface::hasToEvaluateSurrogateModel() const
{
  return 0;
}

bool SolverInterface::hasToEvaluateFineModel() const
{
  return 0;
}

int SolverInterface:: setMeshVertex
(
  int           meshID,
  const double* position )
{
  return -1;
}

int SolverInterface:: getMeshVertexSize
(
  int meshID) const
{
  return -1;
}

void SolverInterface:: setMeshVertices
(
  int           meshID,
  int           size,
  const double* positions,
  int*          ids )
{}

void SolverInterface:: getMeshVertices
(
  int        meshID,
  int        size,
  const int* ids,
  double*    positions ) const
{}

void SolverInterface:: getMeshVertexIDsFromPositions
(
  int           meshID,
  int           size,
  const double* positions,
  int*          ids ) const
{}

int SolverInterface:: setMeshEdge
(
  int meshID,
  int firstVertexID,
  int secondVertexID )
{
  return -1;
}

void SolverInterface:: setMeshTriangle
(
  int meshID,
  int firstEdgeID,
  int secondEdgeID,
  int thirdEdgeID )
{}

void SolverInterface:: setMeshTriangleWithEdges
(
  int meshID,
  int firstVertexID,
  int secondVertexID,
  int thirdVertexID )
{}

void SolverInterface:: setMeshQuad
(
  int meshID,
  int firstEdgeID,
  int secondEdgeID,
  int thirdEdgeID,
  int fourthEdgeID )
{}

void SolverInterface:: setMeshQuadWithEdges
(
  int meshID,
  int firstVertexID,
  int secondVertexID,
  int thirdVertexID,
  int fourthVertexID )
{}

void SolverInterface:: mapReadDataTo
(
  int toMeshID )
{}

void SolverInterface:: mapWriteDataFrom
(
  int fromMeshID )
{}


void SolverInterface:: writeBlockVectorData
(
  int     dataID,
  int     size,
  const int*    valueIndices,
  const double* values )
{
  fake_read_write_buffer.clear();
  for(int i = 0; i < size * this->getDimensions(); i++){
      fake_read_write_buffer.push_back(values[i]); 
    }
}

void SolverInterface:: writeVectorData
(
  int           dataID,
  int           valueIndex,
  const double* value )
{
  fake_read_write_buffer.clear();
  for(int i = 0; i < this->getDimensions(); i++){
      fake_read_write_buffer.push_back(value[i]); 
    }
}

void SolverInterface:: writeBlockScalarData
(
  int           dataID,
  int           size,
  const int*    valueIndices,
  const double* values )
{
  fake_read_write_buffer.clear();
  for(int i = 0; i < size; i++){
      fake_read_write_buffer.push_back(values[i]); 
    }
}

void SolverInterface:: writeScalarData
(
  int    dataID,
  int    valueIndex,
  double value )
{
    fake_read_write_buffer.clear();
    fake_read_write_buffer.push_back(value); 
}

void SolverInterface:: readBlockVectorData
(
  int        dataID,
  int        size,
  const int* valueIndices,
  double*    values ) const
{
  for(int i = 0; i < size * this->getDimensions(); i++){
      values[i] = fake_read_write_buffer[i];
    }
}

void SolverInterface:: readVectorData
(
  int     dataID,
  int     valueIndex,
  double* value ) const
{
  for(int i = 0; i < this->getDimensions(); i++){
      value[i] = fake_read_write_buffer[i];
    }
}

void SolverInterface:: readBlockScalarData
(
  int        dataID,
  int        size,
  const int* valueIndices,
  double*    values ) const
{
  for(int i = 0; i<size; i++){
      values[i] = fake_read_write_buffer[i];
    }
}

void SolverInterface:: readScalarData
(
  int     dataID,
  int     valueIndex,
  double& value ) const
{
    value = fake_read_write_buffer[0]; 
}

namespace constants {

const std::string& actionWriteInitialData()
{
    static std::string dummy ("dummy");
    return dummy;
}

const std::string& actionWriteIterationCheckpoint()
{
    static std::string dummy ("dummy");
    return dummy;
}

const std::string& actionReadIterationCheckpoint()
{
    static std::string dummy ("dummy");
    return dummy;
}

} // namespace precice, constants

} // namespace precice
