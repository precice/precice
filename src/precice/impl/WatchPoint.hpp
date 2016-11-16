#pragma once

#include "SharedPointer.hpp"
#include "io/TXTTableWriter.hpp"
#include "mesh/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include <string>
#include <vector>

namespace precice {
  namespace mesh {
    class Vertex;
  }
}

// ----------------------------------------------------------- CLASS_DEFINITION

namespace precice {
namespace impl {

/**
 * @brief Observes and exports coordinates of a point on the geometry.
 */
class WatchPoint
{
public:

   /**
    * @brief Constructor.
    *
    * @param[in] meshToWatch Mesh to be watched, can be empty on construction.
    */
   WatchPoint (
     const Eigen::VectorXd& pointCoords,
     const mesh::PtrMesh&    meshToWatch,
     const std::string&      exportFilename );

   const mesh::PtrMesh& mesh() const;

   const std::string& filename() const;

   /**
    * @brief Initializes the watch point for exporting point data.
    */
   void initialize();

   /**
    * @brief Writes one line with data of the watchpoint into the output file.
    */
   void exportPointData(double time);

private:

   static logging::Logger _log;

   Eigen::VectorXd _point;

   mesh::PtrMesh _mesh;

   io::TXTTableWriter _txtWriter;

   double _shortestDistance;

   std::vector<double> _weights;

   std::vector<mesh::Vertex*> _vertices;

   std::vector<mesh::PtrData> _dataToExport;

   // @bried Holds the information if this processor is the closest
   bool _isClosest;

   void getValue (
     Eigen::VectorXd&  value,
     mesh::PtrData&    data );

   void getValue (
     double&        value,
     mesh::PtrData& data );
};

}} // namespace precice, impl
