#include "TransformToVRML.hpp"

#include <iostream>
#include <fstream>
#include <vector>

#include "IGES2Vrml.hpp"
#include "STEP2Vrml.hpp"

#include  "TransformToVRML.hpp"

int main ( int argc, char** argv )
{
  if ( argc != 3 ) {
    std::cout << " Usage: ./DataConverter inputFileName outputFilename" << std::endl;
    abort ();
  }

  std::string inputFileName ( argv[1] );
  std::string outputFileName ( argv[2] );

  TransformToVRML myTransformer( inputFileName, outputFileName );

  if ( myTransformer.parseFileType () ) {
    if ( myTransformer.doTransform() ) {
      std::cout << " Transform to VRML succeed!" << std::endl;
    }
    else {
      std::cout << " Transform to VRML failed!" << std::endl;
    }
  }
  else {
    std::cout << "Input file type must be one of: .stp .step .igs .iges" << std::endl;
    return 0;
  }
  return 1;
}


