// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_PLOTTER_GRIDDATA_UNSTRUCTURED_UNSTRUCTURED_GRID_WRITER_H_
#define _TARCH_PLOTTER_GRIDDATA_UNSTRUCTURED_UNSTRUCTURED_GRID_WRITER_H_


#include "tarch/plotter/griddata/Writer.h"


namespace tarch {
  namespace plotter {
    namespace griddata {
      namespace unstructured {
        class UnstructuredGridWriter;
      }
    }
  }
}


class tarch::plotter::griddata::unstructured::UnstructuredGridWriter:
  public tarch::plotter::griddata::Writer {
  public:
    /**
     * This is the vertex writer you have to create to plot the vertices.
     * Besides the pure syntax management, the writer also provides a number
     * generator which provides you with a unique id for each vertex.
     *
     * Please ensure that you call close() on the vertex writer before you
     * close the underlying VTKTextFileWriter.
     */
    class VertexWriter {
      public:
        virtual ~VertexWriter() {}

        virtual int plotVertex(const tarch::la::Vector<2,double>& position) = 0;
        virtual int plotVertex(const tarch::la::Vector<3,double>& position) = 0;

        virtual void close() = 0;
    };

    /**
     * Writes the cell.
     */
    class CellWriter {
      public:
        virtual ~CellWriter() {}

        virtual int plotHexahedron(int vertexIndex[8]) = 0;
        virtual int plotQuadrangle(int vertexIndex[4]) = 0;
        virtual int plotLine(int vertexIndex[2]) = 0;
        virtual int plotPoint(int vertexIndex) = 0;
        virtual int plotTriangle(int vertexIndex[3]) = 0;

        virtual void close() = 0;
    };

    /**
     * Caller has to destroy this instance manually. Do not create more than
     * one vertex writer.
     */
    virtual VertexWriter*   createVertexWriter() = 0;

    /**
     * Caller has to destroy this instance manually. Do not create more than one
     * cell writer.
     */
    virtual CellWriter*     createCellWriter() = 0;
};


#endif
