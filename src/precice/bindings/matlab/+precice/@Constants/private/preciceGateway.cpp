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
                result = actionWriteInitialData();
                break;
            case 1:
                result = actionWriteIterationCheckpoint();
                break;
            case 2:
                result = actionReadIterationCheckpoint();
                break;
            default:
                myMexPrint("An unknown ID was passed.");
                return;
        }
        
        outputs[0] = factory.createCharArray(result);
    }
};