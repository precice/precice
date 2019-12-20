#ifndef PRECICE_DATACONVERTER_IGES2VRML_H_
#define PRECICE_DATACONVERTER_IGES2VRML_H_

#include <IGESControl_Reader.hxx> // Import -Tools
#include <Standard_TypeDef.hxx>
#include <TColStd_HSequenceOfTransient.hxx>
#include <TopTools_HSequenceOfShape.hxx>
#include <TopoDS_Shape.hxx>
#include <VrmlAPI_Writer.hxx> // export -Tools
#include <config.h>           // occ variable for precompiler
#include <string>

/**
 * @brief Interface for visitors, performing transform of CAD files.
 *
 * Translation from IGES to VRML (shaded representation)
 *
 * This interface needs libraries TKIGES, TKVRML of OPENCASCADE,
 * and dl library of C/C++ Language.
 *
 * @author Yuan Fei, Beteuer: Bernhard Gatzhammer
 */
class IGES2Vrml {
public:
  /**
    * @brief Constructor
    */
  IGES2Vrml(void);

  /**
    * @brief Destructor
    */
  virtual ~IGES2Vrml(void);

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
  IGESControl_Reader myReader;

  // @brief path and name of input file, must be set before transform
  Standard_CString myIGESFileName;
};

#endif /* IGES2VRML_H_ */
