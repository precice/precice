// Gateway MexFunction object
#include "mex.hpp"
#include "mexAdapter.hpp"
#include <iostream>
#include <sstream>
#include "precice/SolverInterface.hpp"

using namespace matlab::data;
using matlab::mex::ArgumentList;
using namespace precice;
using namespace precice::constants;

enum class FunctionID {
    _constructor_ = 0,
    _destructor_ = 1,
    configure = 2,
    
    initialize = 10,
    initializeData = 11,
    advance = 12,
    finalize = 13,
    
    getDimensions = 20,
    isCouplingOngoing = 21,
    isReadDataAvailable = 22,
    isWriteDataRequired = 23,
    isTimestepComplete = 24,
    hasToEvaluateSurrogateModel = 25,
    hasToEvaluateFineModel = 26,
    
    isActionRequired = 30,
    fulfilledAction = 31,
    
    hasMesh = 40,
    getMeshID = 41,
    getMeshIDs = 42,
    getMeshHandle = 43,
    setMeshVertex = 44,
    getMeshVertexSize = 45,
    setMeshVertices = 46,
    getMeshVertices = 47,
    getMeshVertexIDsFromPositions = 48,
    setMeshEdge = 49,
    setMeshTriangle = 50,
    setMeshTriangleWithEdges = 51,
    setMeshQuad = 52,
    setMeshQuadWithEdges = 53,
    
    hasData = 60,
    getDataID = 61,
    mapReadDataTo = 62,
    mapWriteDataFrom = 63,
    writeBlockVectorData = 64,
    writeVectorData = 65,
    writeBlockScalarData = 66,
    writeScalarData = 67,
    readBlockVectorData = 68,
    readVectorData = 69,
    readBlockScalarData = 70,
    readScalarData = 71
};

class MexFunction: public matlab::mex::Function {
private:
    SolverInterface* interface;
    ArrayFactory factory;
    bool constructed;
    std::shared_ptr<matlab::engine::MATLABEngine> matlabPtr = getEngine();
    
    void myMexPrint(const std::string text) {
        matlabPtr->feval(u"fprintf",0,std::vector<Array>({factory.createScalar(text)}));
    }
    
public:
    MexFunction(): constructed{false}, factory{}, interface{NULL} {}

    void operator()(ArgumentList outputs, ArgumentList inputs) {
        // Get the function ID from the input
        TypedArray<uint8_t> functionIDArray = inputs[0];
        FunctionID functionID = static_cast<FunctionID>(static_cast<int>(functionIDArray[0]));
        
        // Abort if constructor was not called before, or if constructor 
        // was called on an existing solverInterface
        if (functionID==FunctionID::_constructor_ && constructed) {
            myMexPrint("Constructor was called but interface is alread construced.");
            return;
        }
        if (!constructed && functionID!=FunctionID::_constructor_) {
            myMexPrint("Interface was not construced before.");
            return;
        }
        
        switch (functionID) {
            // 0-9: Construction and Configuration
            // 10-19: Steering
            // 20-29: Status Queries
            // 30-39: Action Methods
            // 40-59: Mesh Access
            // 60-79: Data Access
            case FunctionID::_constructor_:
            {
                const StringArray solverName = inputs[1];
                interface = new SolverInterface(solverName[0],0,1);
                constructed = true;
                break;
            }
            case FunctionID::_destructor_:
            {
                delete interface;
                constructed = false;
                break;
            }
            case FunctionID::configure:
            {
                const StringArray configFileName = inputs[1];
                interface->configure(configFileName[0]);
                break;
            }
            
            case FunctionID::initialize:
            {
                double dt = interface->initialize();
                outputs[0] = factory.createArray<double>({1,1}, {dt});
                break;
            }
            case FunctionID::initializeData: //initializeData
            {
                interface->initializeData();
                break;
            }
            case FunctionID::advance:
            {
                const TypedArray<double> dt_old = inputs[1];
                double dt = interface->advance(dt_old[0]);
                outputs[0] = factory.createArray<double>({1,1}, {dt});
                break;
            }
            case FunctionID::finalize:
            {
                interface->finalize();
                break;
            }
            
            case FunctionID::getDimensions:
            {
                int dims = interface->getDimensions();
                outputs[0] = factory.createArray<uint8_t>({1,1}, {(uint8_t) dims});
                break;
            }
            case FunctionID::isCouplingOngoing:
            {
                bool result = interface->isCouplingOngoing();
                outputs[0] = factory.createArray<bool>({1,1}, {result});
                break;
            }
            case FunctionID::isReadDataAvailable:
            {
                bool result = interface->isReadDataAvailable();
                outputs[0] = factory.createArray<bool>({1,1}, {result});
                break;
            }
            case FunctionID::isWriteDataRequired:
            {
                const TypedArray<double> dt_old = inputs[1];
                bool result = interface->isWriteDataRequired(dt_old[0]);
                outputs[0] = factory.createArray<bool>({1,1}, {result});
                break;
            }
            case FunctionID::isTimestepComplete:
            {
                bool result = interface->isTimestepComplete();
                outputs[0] = factory.createArray<bool>({1,1}, {result});
                break;
            }
            case FunctionID::hasToEvaluateSurrogateModel:
            {
                bool result = interface->hasToEvaluateSurrogateModel();
                outputs[0] = factory.createArray<bool>({1,1}, {result});
                break;
            }
            case FunctionID::hasToEvaluateFineModel:
            {
                bool result = interface->hasToEvaluateFineModel();
                outputs[0] = factory.createArray<bool>({1,1}, {result});
                break;
            }
            
            case FunctionID::isActionRequired:
            {
                const StringArray action = inputs[1];
                bool result = interface->isActionRequired(action[0]);
                outputs[0] = factory.createArray<bool>({1,1}, {result});
                break;
            }
            case FunctionID::fulfilledAction:
            {
                const StringArray action = inputs[1];
                interface->fulfilledAction(action[0]);
                break;
            }
            
            case FunctionID::hasMesh:
            {
                const StringArray meshName = inputs[1];
                bool output = interface->hasMesh(meshName[0]);
                outputs[0] = factory.createScalar<bool>(output);
                break;
            }
            case FunctionID::getMeshID:
            {
                const StringArray meshName = inputs[1];
                int id = interface->getMeshID(meshName[0]);
                outputs[0] = factory.createScalar<int32_t>(id);
                break;
            }
            case FunctionID::getMeshIDs:
            {
                const std::set<int> ids = interface->getMeshIDs();
                outputs[0] = factory.createArray<int32_t>({1,ids.size()}, &*ids.begin(), &*ids.end());
                break;
            }
            
            //getMeshHandle: Not implemented yet
            
            case FunctionID::setMeshVertex:
            {
                const TypedArray<int32_t> meshID = inputs[1];
                const TypedArray<double> position = inputs[2];
                int id = interface->setMeshVertex(meshID[0],&*position.begin());
                outputs[0] = factory.createScalar<int32_t>(id);
                break;
            }
            case FunctionID::getMeshVertexSize:
            {
                const TypedArray<int32_t> meshID = inputs[1];
                int size = interface->getMeshVertexSize(meshID[0]);
                outputs[0] = factory.createScalar<int32_t>(size);
                break;
            }
            case FunctionID::setMeshVertices:
            {
                const TypedArray<int32_t> meshID = inputs[1];
                const TypedArray<int32_t> size = inputs[2];
                const TypedArray<double> positions = inputs[3];
                buffer_ptr_t<int32_t> ids_ptr = factory.createBuffer<int32_t>(size[0]);
                int32_t* ids = ids_ptr.get();
                interface->setMeshVertices(meshID[0],size[0],&*positions.begin(),ids);
                outputs[0] = factory.createArrayFromBuffer<int32_t>({1,size[0]}, std::move(ids_ptr));
                break;
            }
            case FunctionID::getMeshVertices:
            {
                const TypedArray<int32_t> meshID = inputs[1];
                const TypedArray<int32_t> size = inputs[2];
                const TypedArray<int32_t> ids = inputs[3];
                buffer_ptr_t<double> positions_ptr = factory.createBuffer<double>(size[0]);
                double* positions = positions_ptr.get();
                interface->getMeshVertices(meshID[0],size[0],&*ids.begin(),positions);
                outputs[0] = factory.createArrayFromBuffer<double>({1,size[0]}, std::move(positions_ptr));
                break;
            }
            case FunctionID::getMeshVertexIDsFromPositions:
            {
                const TypedArray<int32_t> meshID = inputs[1];
                const TypedArray<int32_t> size = inputs[2];
                const TypedArray<double> positions = inputs[3];
                buffer_ptr_t<int32_t> ids_ptr = factory.createBuffer<int32_t>(size[0]);
                int32_t* ids = ids_ptr.get();
                interface->getMeshVertexIDsFromPositions(meshID[0],size[0],&*positions.begin(),ids);
                outputs[0] = factory.createArrayFromBuffer<int32_t>({1,size[0]}, std::move(ids_ptr));
                break;
            }
            case FunctionID::setMeshEdge:
            {
                const TypedArray<int32_t> meshID = inputs[1];
                const TypedArray<int32_t> firstVertexID = inputs[2];
                const TypedArray<int32_t> secondVertexID = inputs[3];
                int id = interface->setMeshEdge(meshID[0],firstVertexID[0],secondVertexID[0]);
                outputs[0] = factory.createScalar<int32_t>(id);
                break;
            }
            case FunctionID::setMeshTriangle:
            {
                const TypedArray<int32_t> meshID = inputs[1];
                const TypedArray<int32_t> firstEdgeID = inputs[2];
                const TypedArray<int32_t> secondEdgeID = inputs[3];
                const TypedArray<int32_t> thirdEdgeID = inputs[3];
                interface->setMeshTriangle(meshID[0],firstEdgeID[0],secondEdgeID[0],thirdEdgeID[0]);
                break;
            }
            case FunctionID::setMeshTriangleWithEdges:
            {
                const TypedArray<int32_t> meshID = inputs[1];
                const TypedArray<int32_t> firstVertexID = inputs[2];
                const TypedArray<int32_t> secondVertexID = inputs[3];
                const TypedArray<int32_t> thirdVertexID = inputs[3];
                interface->setMeshTriangleWithEdges(meshID[0],firstVertexID[0],secondVertexID[0],thirdVertexID[0]);
                break;
            }
            case FunctionID::setMeshQuad:
            {
                const TypedArray<int32_t> meshID = inputs[1];
                const TypedArray<int32_t> firstEdgeID = inputs[2];
                const TypedArray<int32_t> secondEdgeID = inputs[3];
                const TypedArray<int32_t> thirdEdgeID = inputs[3];
                const TypedArray<int32_t> fourthEdgeID = inputs[4];
                interface->setMeshQuad(meshID[0],firstEdgeID[0],secondEdgeID[0],thirdEdgeID[0],fourthEdgeID[0]);
                break;
            }
            case FunctionID::setMeshQuadWithEdges:
            {
                const TypedArray<int32_t> meshID = inputs[1];
                const TypedArray<int32_t> firstVertexID = inputs[2];
                const TypedArray<int32_t> secondVertexID = inputs[3];
                const TypedArray<int32_t> thirdVertexID = inputs[3];
                const TypedArray<int32_t> fourthVertexID = inputs[3];
                interface->setMeshQuadWithEdges(meshID[0],firstVertexID[0],secondVertexID[0],thirdVertexID[0],fourthVertexID[0]);
                break;
            }
            
            case FunctionID::hasData:
            {
                const StringArray dataName = inputs[1];
                const TypedArray<int32_t> meshID = inputs[2];
                bool output = interface->hasData(dataName[0],meshID[0]);
                outputs[0] = factory.createScalar<bool>(output);
                break;
            }
            case FunctionID::getDataID:
            {
                const StringArray dataName = inputs[1];
                const TypedArray<int32_t> meshID = inputs[2];
                int id = interface->getDataID(dataName[0],meshID[0]);
                outputs[0] = factory.createScalar<int32_t>(id);
                break;
            }
            case FunctionID::mapReadDataTo:
            {
                const TypedArray<int32_t> toMeshID = inputs[1];
                interface->mapReadDataTo(toMeshID[0]);
                break;
            }
            case FunctionID::mapWriteDataFrom:
            {
                const TypedArray<int32_t> fromMeshID = inputs[1];
                interface->mapWriteDataFrom(fromMeshID[0]);
                break;
            }
            case FunctionID::writeBlockVectorData:
            {
                const TypedArray<int32_t> dataID = inputs[1];
                const TypedArray<int32_t> size = inputs[2];
                const TypedArray<int32_t> vertexIDs = inputs[3];
                const TypedArray<double> values = inputs[4];
                interface->writeBlockVectorData(dataID[0],size[0],&*vertexIDs.begin(),&*values.begin());
                break;
            }
            case FunctionID::writeVectorData:
            {
                const TypedArray<int32_t> dataID = inputs[1];
                const TypedArray<int32_t> valueIndex = inputs[2];
                const TypedArray<double> value = inputs[3];
                interface->writeVectorData(dataID[0],valueIndex[0],&*value.begin());
                break;
            }
            case FunctionID::writeBlockScalarData:
            {
                const TypedArray<int32_t> dataID = inputs[1];
                const TypedArray<int32_t> size = inputs[2];
                const TypedArray<int32_t> vertexIDs = inputs[3];
                const TypedArray<double> values = inputs[4];
                interface->writeBlockScalarData(dataID[0],size[0],&*vertexIDs.begin(),&*values.begin());
                break;
            }
            case FunctionID::writeScalarData:
            {
                const TypedArray<int32_t> dataID = inputs[1];
                const TypedArray<int32_t> valueIndex = inputs[2];
                const TypedArray<double> value = inputs[3];
                interface->writeScalarData(dataID[0],valueIndex[0],value[0]);
                break;
            }
            case FunctionID::readBlockVectorData:
            {
                const TypedArray<int32_t> dataID = inputs[1];
                const TypedArray<int32_t> size = inputs[2];
                const TypedArray<int32_t> vertexIDs = inputs[3];
                int32_t dim = interface->getDimensions();
                buffer_ptr_t<double> values_ptr = factory.createBuffer<double>(size[0]*dim);
                double* values = values_ptr.get();
                interface->readBlockVectorData(dataID[0],size[0],&*vertexIDs.begin(),values);
                outputs[0] = factory.createArrayFromBuffer<double>({dim,size[0]}, std::move(values_ptr));
                break;
            }
            case FunctionID::readVectorData:
            {
                const TypedArray<int32_t> dataID = inputs[1];
                const TypedArray<int32_t> valueIndex = inputs[2];
                int32_t dim = interface->getDimensions();
                buffer_ptr_t<double> value_ptr = factory.createBuffer<double>(dim);
                double* value = value_ptr.get();
                interface->readVectorData(dataID[0],valueIndex[0],value);
                outputs[0] = factory.createArrayFromBuffer<double>({dim,1}, std::move(value_ptr));
                break;
            }
            case FunctionID::readBlockScalarData:
            {
                const TypedArray<int32_t> dataID = std::move(inputs[1]);
                const TypedArray<int32_t> size = std::move(inputs[2]);
                const TypedArray<int32_t> valueIndices = std::move(inputs[3]);
                const TypedArray<bool> transpose = std::move(inputs[4]);
                int32_t sizeA, sizeB;
                if (transpose[0]) {
                    sizeA = size[0];
                    sizeB = 1;
                }
                else {
                    sizeA = 1;
                    sizeB = size[0];
                }
                buffer_ptr_t<double> values_ptr = factory.createBuffer<double>(size[0]);
                double* values = values_ptr.get();
                interface->readBlockScalarData(dataID[0],size[0],&*valueIndices.begin(),values);
                outputs[0] = factory.createArrayFromBuffer<double>({sizeA,sizeB}, std::move(values_ptr));
                break;
            }
            case FunctionID::readScalarData:
            {
                const TypedArray<int32_t> dataID = inputs[1];
                const TypedArray<int32_t> valueIndex = inputs[2];
                double value;
                interface->readScalarData(dataID[0],valueIndex[0],value);
                outputs[0] = factory.createScalar<double>(value);
                break;
            }
            default:
                myMexPrint("An unknown ID was passed.");
                return;
        }
        // Do error handling
        // myMexPrint("A problem occurred while executing the function.");
    }
};
