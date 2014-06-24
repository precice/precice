/*
 * CCAGridArrayWriter.h
 *
 *  Created on: Dec 1, 2010
 *      Author: atanasoa
 */

#ifndef CCAGRIDARRAYWRITER_H_
#define CCAGRIDARRAYWRITER_H_

#include "tarch/plotter/griddata/regular/CartesianGridArrayWriter.h"


namespace tarch {
  namespace plotter {
    namespace griddata {
      namespace regular {
        namespace cca {
          class CCAGridArrayWriter;
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
class tarch::plotter::griddata::regular::cca::CCAGridArrayWriter:
  public tarch::plotter::griddata::regular::CartesianGridArrayWriter {

  public:
    CCAGridArrayWriter(
      const tarch::la::Vector<2,int>&     numberOfGridPoints,
      const tarch::la::Vector<2,double>&  domainSize,
      const tarch::la::Vector<2,double>&  origin
    );

    CCAGridArrayWriter(
      const tarch::la::Vector<3,int>&     numberOfGridPoints,
      const tarch::la::Vector<3,double>&  domainSize,
      const tarch::la::Vector<3,double>&  origin
    );

    virtual ~CCAGridArrayWriter();

    /**
     * t.b.d.
     */
    virtual void writeToVertexArray( double* array,const int length);
    virtual void writeToCellArray( double* array,const int length);
    virtual void writeToFile( const std::string& filename );


};


#endif /* CCAGRIDARRAYWRITER_H_ */
