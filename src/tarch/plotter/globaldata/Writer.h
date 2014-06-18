// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_PLOTTER_GLOBAL_DATA_WRITER_H_
#define _TARCH_PLOTTER_GLOBAL_DATA_WRITER_H_

#include <string>
#include "tarch/la/Vector.h"


namespace tarch {
  namespace plotter {
    namespace globaldata {
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
class tarch::plotter::globaldata::Writer{
  public :
    /**
     * A writer for scalar data on elements.
     */
    class DataWriter {
      public:
        virtual ~DataWriter(){}

        virtual void close() = 0;

        /**
         * Write data for one cell.
         *
         * @param index Index of the vertex. This index has to equal the index
         *              used for the cell within the VTKWriter class
         *              interface.
         * @param value Value for the cell.
         */
        virtual void plotRecord( double value ) = 0;
        virtual void plotRecord( const tarch::la::Vector<2,double>& value ) = 0;
        virtual void plotRecord( const tarch::la::Vector<3,double>& value ) = 0;
        virtual void plotRecord( int numberOfEntries, double entries[] ) = 0;

        virtual double getMinValue() const = 0;
        virtual double getMaxValue() const = 0;
        virtual int    getMaxNumberOfEntriesPerRecord() const = 0;
    };

    virtual ~Writer(){}

    /**
     * Close the writer. If you do so, all derived writers such as vertex data
     * writers or cell writers have to be closed.
     */
    virtual void close() = 0;

    /**
     * @return Whether writer is ready to accept data.
     */
    virtual bool isOpen() = 0;

    /**
     * You are responsible to destroy the instance on the heap explicitly.
     */
    virtual DataWriter* createDataWriter( const std::string& dataIdentifier ) = 0;
};

#endif
