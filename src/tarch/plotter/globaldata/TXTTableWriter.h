// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_PLOTTER_GLOBAL_DATA_TXT_TABLE_WRITER_H_
#define _TARCH_PLOTTER_GLOBAL_DATA_TXT_TABLE_WRITER_H_

#include "tarch/plotter/globaldata/Writer.h"

#include <sstream>
#include <fstream>

namespace tarch {
  namespace plotter {
    namespace globaldata {
      class TXTTableWriter;
    }
  }
}


/**
 * Writes data into a text file that is organised as a table.
 *
 * @author Tobias Weinzierl
 */
class tarch::plotter::globaldata::TXTTableWriter: public tarch::plotter::globaldata::Writer {
  private:
    /**
     * Header statement (first line of file).
     */
    static const std::string HEADER;

    /**
     * Header statement (first line of file).
     */
    static const std::string SEPARATOR;

    /**
     * Output stream.
     */
    std::ofstream _out;
  public:
    TXTTableWriter( const std::string& filename );

    virtual ~TXTTableWriter();

    /**
     * Close the writer. If you do so, all derived writers such as vertex data
     * writers or cell writers have to be closed.
     */
    virtual void close();

    /**
     * @return Whether writer is ready to accept data.
     */
    virtual bool isOpen();

    /**
     * You are responsible to destroy the instance on the heap explicitly.
     */
    virtual DataWriter*    createDataWriter( const std::string& dataIdentifier );

    /**
     * A writer for scalar data on elements.
     */
    class DataWriter {
      private:
        TXTTableWriter&   _writer;
      protected:
        DataWriter( TXTTableWriter& writer, const std::string& dataIdentifier );

      public:
        virtual ~DataWriter();

        virtual void close() = 0;

        /**
         * Write data for one cell.
         *
         * @param index Index of the vertex. This index has to equal the index
         *              used for the cell within the VTKWriter class
         *              interface.
         * @param value Value for the cell.
         */
        virtual void plotRecord( int index, double value ) = 0;
        virtual void plotRecord( int index, const tarch::la::Vector<2,double>& value ) = 0;
        virtual void plotRecord( int index, const tarch::la::Vector<3,double>& value ) = 0;

        virtual tarch::la::Vector<3,double> getMinValue() const = 0;
        virtual tarch::la::Vector<3,double> getMaxValue() const = 0;
    };
};

#endif
