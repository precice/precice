// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_PLOTTER_GRIDDATA_REGULAR_VTK_TEXTFILEWRITER_H_
#define _TARCH_PLOTTER_GRIDDATA_REGULAR_VTK_TEXTFILEWRITER_H_


#include "tarch/plotter/griddata/regular/CartesianGridArrayWriter.h"


namespace tarch {
  namespace plotter {
    namespace griddata {
      namespace regular {
        namespace vtk {
          class VTKTextFileWriter;
        }
      }
    }
  }
}


/**
 * VTK Writer
 *
 *
 * @author Tobias Weinzierl, Atanas Atanasov
 * @version $Revision: 1.1 $
 */
class tarch::plotter::griddata::regular::vtk::VTKTextFileWriter:
  public tarch::plotter::griddata::regular::CartesianGridArrayWriter {
  private:
    static const std::string HEADER;
    const int _precision;
    const std::string _doubleOrFloat;

    std::string setDoubleOrFloatString(const int precision){
      if (precision < 7){
        return "float";
      } else {
        return "double";
      }
    }


  public:
    VTKTextFileWriter(
      const tarch::la::Vector<2,int>&     numberOfGridPoints,
      const tarch::la::Vector<2,double>&  domainSize,
      const tarch::la::Vector<2,double>&  origin,
      const int precision=6
    );

    VTKTextFileWriter(
      const tarch::la::Vector<3,int>&     numberOfGridPoints,
      const tarch::la::Vector<3,double>&  domainSize,
      const tarch::la::Vector<3,double>&  origin,
      const int precision=6
    );

    virtual ~VTKTextFileWriter();

    /**
     * t.b.d.
     */
    virtual void writeToFile( const std::string& filename );
};

#endif
