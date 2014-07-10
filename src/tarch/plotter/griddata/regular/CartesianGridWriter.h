// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_PLOTTER_GRIDDATA_REGULAR_CARTESIAN_GRID_WRITER_H_
#define _TARCH_PLOTTER_GRIDDATA_REGULAR_CARTESIAN_GRID_WRITER_H_


#include "tarch/plotter/griddata/Writer.h"


namespace tarch {
  namespace plotter {
    namespace griddata {
      namespace regular {
        class CartesianGridWriter;
      }
    }
  }
}


/**
 * Cartesian Grid Writer
 *
 * Abstract Interface for all the Cartesian grid writers. Provides access
 * operations for a user to find out which vertex is the right one if you
 * wanna plot a vertex's data and if you know only the vertex position.
 *
 * @author Atanas Atanasov, Tobias Weinzierl
 */
class tarch::plotter::griddata::regular::CartesianGridWriter:
  public tarch::plotter::griddata::Writer {
  public:
	/**
	 * Which Vertex is Near to Position?
	 *
	 * This operation searches for the nearest vertex of position and returns
	 * its index. The result is always greater equal zero.
	 */
    virtual int getVertexIndex( const tarch::la::Vector<2,double>& position ) = 0;

	/**
	 * Which Vertex is Near to Position?
	 *
	 * This operation searches for the nearest vertex of position and returns
	 * its index.
	 */
    virtual int getVertexIndex( const tarch::la::Vector<3,double>& position ) = 0;

	/**
	 * Which Vertex is Near to Position?
	 *
	 * This operation searches for the cell of position and returns
	 * its index.
	 *
	 * @param position offset of cell
	 */
    virtual int getCellIndex( const tarch::la::Vector<2,double>& position ) = 0;

	/**
	 * Which Cell is Near to Position?
	 *
	 * This operation searches for the nearest vertex of position and returns
	 * its index. The result is always greater equal zero.
	 */
    virtual int getCellIndex( const tarch::la::Vector<3,double>& position ) = 0;

    /**
     * Does the element specified contain any vertices?
     */
    virtual bool containsVertex(
      const tarch::la::Vector<3,double>& position
    ) const = 0;

    /**
     * Does the element specified contain any vertices?
     */
    virtual bool containsVertex(
      const tarch::la::Vector<2,double>& position
    ) const = 0;

    /**
     * Does the element specified contain any cells?
     */
    virtual bool containsCell(
      const tarch::la::Vector<2,double>& offset,
      const tarch::la::Vector<2,double>& boundingBox
    ) const = 0;

    /**
     * Does the element specified contain any cells?
     */
    virtual bool containsCell(
      const tarch::la::Vector<3,double>& offset,
      const tarch::la::Vector<3,double>& boundingBox
    ) const = 0;

    virtual tarch::la::Vector<3,double> getH() const = 0;

    /**
     * A very simple extension of the writer's interface that provides three
     * proxy operations. Hence, this interface does not introduce new features,
     * it just gives access to the underlying (parent) class' operations
     * through the data writer interface.
     */
    class VertexDataWriter: public tarch::plotter::griddata::Writer::VertexDataWriter {
      public:
        virtual tarch::la::Vector<3,double> getH() const = 0;
        virtual int getVertexIndex( const tarch::la::Vector<2,double>& position ) = 0;
        virtual int getVertexIndex( const tarch::la::Vector<3,double>& position ) = 0;
        virtual bool containsVertex( const tarch::la::Vector<3,double>& position ) const = 0;
        virtual bool containsVertex( const tarch::la::Vector<2,double>& position ) const = 0;
    };

    /**
     * A very simple extension of the writer's interface that provides three
     * proxy operations. Hence, this interface does not introduce new features,
     * it just gives access to the underlying (parent) class' operations
     * through the data writer interface.
     */
    class CellDataWriter: public tarch::plotter::griddata::Writer::CellDataWriter {
      public:
        virtual tarch::la::Vector<3,double> getH() const = 0;
        virtual int getCellIndex( const tarch::la::Vector<2,double>& position ) = 0;
        virtual int getCellIndex( const tarch::la::Vector<3,double>& position ) = 0;
        virtual bool containsCell( const tarch::la::Vector<2,double>& offset, const tarch::la::Vector<2,double>& boundingBox ) const = 0;
        virtual bool containsCell( const tarch::la::Vector<3,double>& offset, const tarch::la::Vector<3,double>& boundingBox ) const = 0;
    };
};


#endif
