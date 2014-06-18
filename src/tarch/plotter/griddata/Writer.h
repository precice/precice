// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_PLOTTER_GRIDDATA_WRITER_H_
#define _TARCH_PLOTTER_GRIDDATA_WRITER_H_


#include "tarch/la/Vector.h"


namespace tarch {
  namespace plotter {
    namespace griddata {
  	  class Writer;
    }
  }
}


/**
 * Interface for all Plotters.
 *
 * @author Atanas Atanasov
 * @author Tobias Neckel
 * @author Tobias Weinzierl
 */
class tarch::plotter::griddata::Writer {
  private:
    /**
     * Copy constructor
     *
     * A writer never should be copied. However, we have to make this operation
     * protected to allow implementations to hide their copy constructor as
     * well.
     */
    Writer(const Writer& writer){}

    /**
     * Assignment operator.
     *
     * A writer never should be copied. However, we have to make this operation
     * protected to allow implementations to hide their copy constructor as
     * well.
     */
    Writer& operator=(const Writer& writer) {return *this;}

	public:
    Writer() {}

    virtual ~Writer() {}

    virtual void writeToFile( const std::string& filename ) = 0;

    /**
     * @return Whether writer is ready to accept data.
     */
    virtual bool isOpen() = 0;

    /**
     * Clear the writer, i.e. erase all the data. However, as the writer does
     * not track how many vertex and cell writers you've created, it's up to
     * you to ensure that none of these instances is left.
     */
    virtual void clear() = 0;

    /**
     * A writer for scalar data on elements.
     */
    class CellDataWriter {
   	  public:
        virtual ~CellDataWriter() {};

        /**
         * Write data for one cell.
         *
         * @param index Index of the cell. This index has to equal the index
         *              used for the cell within the VTKWriter class
         *              interface.
         * @param value Value for the cell.
         */
        virtual void plotCell( int index, double value ) = 0;
        virtual void plotCell( int index, const tarch::la::Vector<2,double>& value ) = 0;
        virtual void plotCell( int index, const tarch::la::Vector<3,double>& value ) = 0;

        virtual double getMinValue() const = 0;
        virtual double getMaxValue() const = 0;

        virtual void close() = 0;
    };

    /**
     * A writer for scalar data on points (vertices).
     */
    class VertexDataWriter {
      public:
        virtual ~VertexDataWriter() {};

        /**
         * Write data for one cell.
         *
         * @param index Index of the vertex. This index has to equal the index
         *              used for the cell within the VTKWriter class
         *              interface.
         * @param value Value for the cell.
         */
        virtual void plotVertex( int index, double value ) = 0;
        virtual void plotVertex( int index, const tarch::la::Vector<2,double>& value ) = 0;
        virtual void plotVertex( int index, const tarch::la::Vector<3,double>& value ) = 0;

        virtual double getMinValue() const = 0;
        virtual double getMaxValue() const = 0;

        virtual void close() = 0;
    };

    /**
     * Caller has to destroy this instance manually.
     */
    virtual CellDataWriter*    createCellDataWriter( const std::string& identifier, int recordsPerCell )   = 0;

    /**
     * Caller has to destroy this instance manually.
     */
    virtual VertexDataWriter*  createVertexDataWriter( const std::string& identifier, int recordsPerVertex ) = 0;
};

#endif
