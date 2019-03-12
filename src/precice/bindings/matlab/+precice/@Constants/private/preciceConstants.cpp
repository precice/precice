// Gateway MexFunction object
#include "mex.hpp"
#include "mexAdapter.hpp"
#include <iostream>
#include "precice/Constants.hpp"

using namespace matlab::data;
using matlab::mex::ArgumentList;
using namespace precice::constants;

class MexFunction: public matlab::mex::Function {
private:
    ArrayFactory factory;
    void myMexPrint(const std::string text) {
        std::cout << "MEX constants gateway: " << text << std::endl;
    }

public:
    MexFunction(): factory{} {}

    void operator()(ArgumentList outputs, ArgumentList inputs) {
        // define the constantID
        TypedArray<uint8_t> functionIDArray = inputs[0];
        int constantID = functionIDArray[0];
        std::cout << "Constant " << (int)constantID << " was called." << std::endl;
        std::string result;
        
        // assign
        switch (constantID) {
            case 0:
                result = nameConfiguration();
                break;
            case 1:
                result = dataDisplacements();
                break;
            case 2:
                result = dataForces();
                break;
            case 3:
                result = dataVelocities();
                break;
            case 4:
                result = actionWriteInitialData();
                break;
            case 5:
                result = actionWriteIterationCheckpoint();
                break;
            case 6:
                result = actionReadIterationCheckpoint();
                break;
            case 7:
                result = actionPlotOutput();
                break;
            case 8:
                result = exportVTK();
                break;
            case 9:
                result = exportAll();
                break;
            default:
                myMexPrint("An unknown ID was passed.");
                return;
        }
        
        outputs[0] = factory.createCharArray(result);
    }
};