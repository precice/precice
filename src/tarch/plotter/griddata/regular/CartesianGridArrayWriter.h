// Copyright (C) 2009 Technische Universitaet Muenchen
// This file is part of the Peano project. For conditions of distribution and
// use, please see the copyright notice at www5.in.tum.de/peano
#ifndef _TARCH_PLOTTER_GRIDDATA_REGULAR_CARTESIAN_GRID_ARRAY_WRITER_H_
#define _TARCH_PLOTTER_GRIDDATA_REGULAR_CARTESIAN_GRID_ARRAY_WRITER_H_


#include "logging/Logger.hpp"
#include "tarch/plotter/griddata/regular/CartesianGridWriter.h"

#include <vector>

namespace tarch {
  namespace plotter {
    namespace griddata {
      namespace regular {
        class CartesianGridArrayWriter;
      }
    }
  }
}


/**
 * Abstract Cartesian Grid Writer
 *
 * This writer buffers all the output data into an array. When you close the
 * writers, these array buffers have to be written to a data file due to
 * writeToFile() which is the only abstract function of this class.
 * writeToFile() also has to take care to format the output in the right
 * syntax.
 *
 * @author Tobias Weinzierl, Atanas Atanasov
 * @version $Revision: 1.1 $
 */
class tarch::plotter::griddata::regular::CartesianGridArrayWriter:
  public tarch::plotter::griddata::regular::CartesianGridWriter {
  protected:

    /**
     * Logging device.
     */
    static logging::Logger _log;


    struct DataSet {
      std::string _identifier;
      int         _recordsPerEntry;
      double*     _data;
      DataSet( const std::string& identifier, int recordsPerEntry, const tarch::la::Vector<3,int>& numberOfGridPoints );
      void clear();
    };

    bool _writtenToFile;

    tarch::la::Vector<3,int>     _numberOfGridPoints;
    tarch::la::Vector<3,double>  _domainSize;
    tarch::la::Vector<3,double>  _origin;

    std::vector<DataSet> _vertexData;
    std::vector<DataSet> _cellData;

    tarch::la::Vector<3,int> getNumberOfCells() const;

  public:
    CartesianGridArrayWriter(
      const tarch::la::Vector<2,int>&     numberOfGridPoints,
      const tarch::la::Vector<2,double>&  domainSize,
      const tarch::la::Vector<2,double>&  origin
    );

    CartesianGridArrayWriter(
      const tarch::la::Vector<3,int>&     numberOfGridPoints,
      const tarch::la::Vector<3,double>&  domainSize,
      const tarch::la::Vector<3,double>&  origin
    );

    virtual ~CartesianGridArrayWriter();

    virtual tarch::la::Vector<3,double> getH() const;

    /**
     * From the interface Writer.
     */
    virtual bool isOpen();

    /**
     * From the interface Writer.
     */
    virtual void clear();

    /**
     * From the interface Writer. The result type actually is a
     * CartesianGridWriter::CellDataWriter, so if you'd like to, you might
     * downcast it.
     */
    virtual tarch::plotter::griddata::Writer::CellDataWriter*    createCellDataWriter( const std::string& identifier, int recordsPerCell );

    /**
     * From the interface Writer. The result type actually is a
     * CartesianGridWriter::VertexDataWriter, so if you'd like to, you might
     * downcast it.
     */
    virtual tarch::plotter::griddata::Writer::VertexDataWriter*  createVertexDataWriter( const std::string& identifier, int recordsPerVertex );

    /**
     * From the CartesianGridWriter.
     */
    virtual int getVertexIndex( const tarch::la::Vector<2,double>& position );

    /**
     * From the CartesianGridWriter.
     *
     * Implementation details:
     * - Move vertex @f$ \frac{h}{2} @f$ to reconstruct the rounding effect.
     * - Then
     */
    virtual int getVertexIndex( const tarch::la::Vector<3,double>& position );

    /**
     * From the CartesianGridWriter.
     */
    virtual int getCellIndex( const tarch::la::Vector<2,double>& position );

    /**
     * From the CartesianGridWriter.
     */
    virtual int getCellIndex( const tarch::la::Vector<3,double>& position );

    /**
     * From the CartesianGridWriter.
     *
     * The operation returns false
     * - if the voxel specified does not contain a cell,
     * - if
     */
    virtual bool containsVertex(
      const tarch::la::Vector<3,double>& x
    ) const;

    /**
     * From the CartesianGridWriter.
     *
     * Converts both 2d-vectors into threedimensional vectors and then forwards
     * them to the other containsVertex() operation.
     */
    virtual bool containsVertex(
      const tarch::la::Vector<2,double>& x
    ) const;

    /**
     * From the CartesianGridWriter.
     *
     * Converts both 2d-vectors into threedimensional vectors and then forwards
     * them to the other containsCell() operation.
     */
    virtual bool containsCell(
      const tarch::la::Vector<2,double>& offset,
      const tarch::la::Vector<2,double>& boundingBox
    ) const;

    /**
     * From the CartesianGridWriter.
     *
     * Basically, this operation evaluates a cut of the voxel passed with the
     * domain. Fot his, we follow the classical Minkovski-sum paradigm:
     * - Compute the center and half the mesh size of the voxel passed.
     * - Extend the computational domain of the plotter by this half mesh size.
     * - Check whether the center of the voxel (point) is inside this extended
     *   domain.
     */
    virtual bool containsCell(
      const tarch::la::Vector<3,double>& offset,
      const tarch::la::Vector<3,double>& boundingBox
    ) const;


    class CellDataWriter: public tarch::plotter::griddata::regular::CartesianGridWriter::CellDataWriter {
      private:
        /**
         * The father class is a friend. There are no other friends.
         */
        friend class CartesianGridArrayWriter;

        /**
         * Underlying VTK writer.
         */
        DataSet& _dataSet;
        CartesianGridArrayWriter& _parent;

        double _minValue;
        double _maxValue;

        CellDataWriter(DataSet& writer, CartesianGridArrayWriter& parent);
      public:
        virtual ~CellDataWriter();

        virtual void close();

        virtual void plotCell( int index, double value );
        virtual void plotCell( int index, const tarch::la::Vector<2,double>& value );
        virtual void plotCell( int index, const tarch::la::Vector<3,double>& value );

        virtual double getMinValue() const;
        virtual double getMaxValue() const;

        virtual tarch::la::Vector<3,double> getH() const;
        virtual int getCellIndex( const tarch::la::Vector<2,double>& position );
        virtual int getCellIndex( const tarch::la::Vector<3,double>& position );
        virtual bool containsCell( const tarch::la::Vector<2,double>& offset, const tarch::la::Vector<2,double>& boundingBox ) const;
        virtual bool containsCell( const tarch::la::Vector<3,double>& offset, const tarch::la::Vector<3,double>& boundingBox ) const;
    };

    class VertexDataWriter: public tarch::plotter::griddata::regular::CartesianGridWriter::VertexDataWriter {
      private:
        /**
         * The father class is a friend. There are no other friends.
         */
        friend class CartesianGridArrayWriter;

        /**
         * Underlying VTK writer.
         */
        DataSet& _dataSet;
        CartesianGridArrayWriter& _parent;

        double _minValue;
        double _maxValue;

        VertexDataWriter(DataSet& writer, CartesianGridArrayWriter& parent);
      public:
        virtual ~VertexDataWriter();

        virtual void close();

        virtual void plotVertex( int index, double value );
        virtual void plotVertex( int index, const tarch::la::Vector<2,double>& value );
        virtual void plotVertex( int index, const tarch::la::Vector<3,double>& value );

        virtual double getMinValue() const;
        virtual double getMaxValue() const;

        virtual tarch::la::Vector<3,double> getH() const;
        virtual int getVertexIndex( const tarch::la::Vector<2,double>& position );
        virtual int getVertexIndex( const tarch::la::Vector<3,double>& position );
        virtual bool containsVertex( const tarch::la::Vector<3,double>& position ) const;
        virtual bool containsVertex( const tarch::la::Vector<2,double>& position ) const;
    };
};

#endif

