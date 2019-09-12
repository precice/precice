// Gateway MexFunction object
#include "mex.hpp"
#include "mexAdapter.hpp"
#include <iostream>
#include "precice/SolverInterface.hpp"

using namespace matlab::data;
using matlab::mex::ArgumentList;
using namespace precice::constants;

class MexFunction: public matlab::mex::Function {
private:
    ArrayFactory factory;

public:
    void operator()(ArgumentList outputs, ArgumentList inputs) {
        // define the constantID
        TypedArray<uint8_t> functionIDArray = inputs[0];
        int constantID = functionIDArray[0];
        std::string result;
        
        // assign
        switch (constantID) {
            case 0:
                result = actionWriteInitialData();
                break;
            case 1:
                result = actionWriteIterationCheckpoint();
                break;
            case 2:
                result = actionReadIterationCheckpoint();
                break;
            default:
                std::cout << "MEX constants gateway: An unknown ID was passed." << std::endl;
                return;
        }
        
        outputs[0] = factory.createCharArray(result);
    }
};
