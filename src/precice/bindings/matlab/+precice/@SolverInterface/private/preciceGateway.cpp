// Gateway MexFunction object
#include "mex.hpp"
#include "mexAdapter.hpp"
#include <iostream>
#include <sstream>
#include "precice/SolverInterface.hpp"
#include "precice/Constants.hpp"

using namespace matlab::data;
using matlab::mex::ArgumentList;
using namespace precice;
using namespace precice::constants;

typedef bool (MexFunction::*fnPtr)(ArgumentList&,ArgumentList&);
struct fnPtrWrap {
    fnPtr functionName;
    bool valid;
};

class MexFunction: public matlab::mex::Function {
private:
    SolverInterface* interface;
    ArrayFactory factory;
    bool constructed;
    
    void myMexPrint(const std::string text) {
        std::cout << "MEX gateway: " << text << std::endl;
    }
    
    bool checkInput(ArgumentList& outputs, const int& noutputs, 
            ArgumentList& inputs, const int& ninputs, const ArrayType* types) {
        if (inputs.size() < ninputs+1 || outputs.size() < noutputs){
            myMexPrint("Wrong number of input or output arguments.");
            return false;
        }
        for (int j = 0; j<ninputs; j++) {
            if (inputs[j+1].getType() != types[j]) {
                myMexPrint("False type of input arguments.");
                return false;
            }
        }
        return true;
    }
    
    bool constructInterface(ArgumentList& inputs, ArgumentList& outputs) {
        CharArray solverName = inputs[1];
        interface = new SolverInterface(solverName.toAscii(),0,1);
        constructed = true;
        return true;
    }
    
    bool destructInterface(ArgumentList& inputs, ArgumentList& outputs) {
        delete interface;
        constructed = false;
        return true;
    }
    
    bool configure(ArgumentList& inputs, ArgumentList& outputs) {
        CharArray configFileName = inputs[1];
        interface->configure(configFileName.toAscii());
        return true;
    }
    
    bool initialize(ArgumentList& inputs, ArgumentList& outputs) {
        double dt = interface->initialize();
        outputs[0] = factory.createArray<double>({1,1}, {dt});
        return true;
    }
    
    bool advance(ArgumentList& inputs, ArgumentList& outputs) {
        TypedArray<double> dt_old = std::move(inputs[1]);
        double dt = interface->advance(dt_old[0]);
        outputs[0] = factory.createArray<double>({1,1}, {dt});
        return true;
    }
    
    bool finalize(ArgumentList& inputs, ArgumentList& outputs) {
        interface->finalize();
        return true;
    }
    
    bool getDimensions(ArgumentList& inputs, ArgumentList& outputs) {
        int dims = interface->getDimensions();
        outputs[0] = factory.createArray<uint8_t>({1,1}, {(uint8_t) dims});
        return true;
    }
    
    bool isCouplingOngoing(ArgumentList& inputs, ArgumentList& outputs) {
        bool result = interface->isCouplingOngoing();
        outputs[0] = factory.createArray<bool>({1,1}, {result});
        return true;
    }
    
    bool isActionRequired(ArgumentList& inputs, ArgumentList& outputs) {
        CharArray action = inputs[1];
        bool result = interface->isActionRequired(action.toAscii());
        outputs[0] = factory.createArray<bool>({1,1}, {result});
        return true;
    }
    
    bool fulfilledAction(ArgumentList& inputs, ArgumentList& outputs) {
        CharArray action = inputs[1];
        interface->fulfilledAction(action.toAscii());
        return true;
    }
    
    bool getMeshID(ArgumentList& inputs, ArgumentList& outputs) {
        CharArray meshName = inputs[1];
        int id = interface->getMeshID(meshName.toAscii());
        outputs[0] = factory.createArray<uint64_t>({1,1}, {(uint64_t)id});
        return true;
    }
    
    bool setMeshVertex(ArgumentList& inputs, ArgumentList& outputs) {
        TypedArray<uint64_t> meshID = std::move(inputs[1]);
        TypedArray<double> position = std::move(inputs[2]);
        int id = interface->setMeshVertex(meshID[0],&*position.begin());
        outputs[0] = factory.createArray<uint64_t>({1,1}, {(uint64_t) id});
        return true;
    }

public:
    MexFunction(): constructed{false}, factory{}, interface{NULL} {}

    void operator()(ArgumentList outputs, ArgumentList inputs) {
        bool valid=false;
        fnPtr functionName=NULL;
        if (inputs.size() == 0){
            myMexPrint("No argument.");
            return;
        }
        if (inputs[0].getType() != ArrayType::UINT8){
            myMexPrint("First argument no function ID.");
            return;
        }
        // define the functionID
        TypedArray<uint8_t> functionIDArray = inputs[0];
        int functionID = functionIDArray[0];
        std::cout << "Gateway function " << (int)functionID << " was called." << std::endl;
        
        // Abort if constructor was not called before, or if constructor 
        // was called on an existing solverInterface
        if (functionID==0 && constructed) {
            myMexPrint("Constructor was called but interface is alread construced.");
            return;
        }
        if (!constructed && functionID!=0) {
            myMexPrint("Interface was not construced before.");
            return;
        }
        
        int ninputs = 0;
        int noutputs = 0;
        ArrayType* types;
        switch (functionID) {
            // 0-9: Construction and Configuration
            // 10-19: Steering
            // 20-29: Status Queries
            // 30-39: Action Methods
            // 40-59: Mesh Access
            // 60-79: Data Access
            case 0: //construction
                ninputs = 1;
                types = new ArrayType [ninputs] {ArrayType::CHAR};
                functionName = &MexFunction::constructInterface;
                break;
            case 1: //destruction
                functionName = &MexFunction::destructInterface;
                break;
            case 2: //configuration
                ninputs = 1;
                types = new ArrayType [ninputs] {ArrayType::CHAR};
                functionName = &MexFunction::configure;
                break;
            case 10: //initialization
                noutputs = 1;
                functionName = &MexFunction::initialize;
                break;
            case 12: //advance
                ninputs = 1;
                types = new ArrayType [ninputs] {ArrayType::DOUBLE};
                noutputs = 1;
                functionName = &MexFunction::advance;
                break;
            case 13: //finalization
                functionName = &MexFunction::finalize;
                break;
            case 20: //getDimensions
                noutputs = 1;
                functionName = &MexFunction::getDimensions;
                break;
            case 21: //isCouplingOngoing
                noutputs = 1;
                functionName = &MexFunction::isCouplingOngoing;
                break;
            case 30: //isActionRequired
                ninputs = 1;
                types = new ArrayType [ninputs] {ArrayType::CHAR};
                noutputs = 1;
                functionName = &MexFunction::isActionRequired;
                break;
            case 31: //fulfilledAction
                ninputs = 1;
                types = new ArrayType [ninputs] {ArrayType::CHAR};
                functionName = &MexFunction::fulfilledAction;
                break;
            case 41: //getMeshID
                ninputs = 1;
                types = new ArrayType [ninputs] {ArrayType::CHAR};
                noutputs = 1;
                functionName = &MexFunction::getMeshID;
                break;
            case 44: //setMeshVertex
                ninputs = 2;
                types = new ArrayType [ninputs] {ArrayType::UINT64,ArrayType::DOUBLE};
                noutputs = 1;
                functionName = &MexFunction::setMeshVertex;
                break;
            default:
                myMexPrint("An unknown ID was passed.");
                return;
        }
        valid = checkInput(outputs,noutputs,inputs,ninputs,types);
        if (valid) {
            bool success = (*this.*functionName)(inputs,outputs);
            if (!success) {
                myMexPrint("A problem occurred while executing the function.");
            }
        }
    }
};