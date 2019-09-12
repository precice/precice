#ifndef PRECICE_DATACONVERTER_STEP2VRML_HPP_
#define PRECICE_DATACONVERTER_STEP2VRML_HPP_

#include <STEPControl_Reader.hxx>           // Import -Tools
#include <Standard_TypeDef.hxx>             //standard data type in Opencascade
#include <TColStd_HSequenceOfTransient.hxx> //save list of shapes for transfer
#include <TopTools_HSequenceOfShape.hxx>
#include <TopoDS_Shape.hxx>
#include <VrmlAPI_Writer.hxx> // export -Tools
#include <config.h>           // occ variable for precompiler
#include <string>

/**
 * @brief Interface for visitors, performing transform of CAD files.
 *
 * Translation from STEP to VRML (shaded representation)
 *
 * This interface needs libraries TKSTEP, TKVRML of OPENCASCADE,
 * and dl library of C/C++ Language.
 *
 * @author Yuan Fei, Beteuer: Bernhard Gatzhammer
 */

class STEP2Vrml {
public:
  /**
    * @brief Constructor
    */
  STEP2Vrml(void);

  /**
    * @brief Destructor
    */
  virtual ~STEP2Vrml(void);

  /**
    * @brief setting input file
    */
  bool readFile(const std::string &inputFileName);

  /**
    * @brief prints Statistics of input file
    */
  void loadReport(void);

  /**
    * @brief transform to VRML
    */
  bool transfer(const std::string &outputFileName);

private:
  // @brief reader for IGES file
  STEPControl_Reader myReader;

  // @brief path and name of input file, must be set before transform
  Standard_CString mySTEPFileName;
};

#endif /* STEP2VRML_H_ */
