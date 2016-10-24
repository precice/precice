#pragma once

#include "spacetree/Spacetree.hpp"
#include <string>
#include <list>

namespace precice {
namespace spacetree {

/**
 * @brief Writes the leave cells of a spacetree as voxels into a vtk file
 */
class ExportSpacetree : public Spacetree::Visitor
{
public:

   /**
    * @brief Constructor, defining vtk filename.
    */
   ExportSpacetree ( const std::string& location, const std::string& filename );

   virtual ~ExportSpacetree() {}

   void doExport ( Spacetree& toExport );

   void nodeCallback (
     const Eigen::VectorXd& center,
     const Eigen::VectorXd& halflengths,
     int                    position );

   void leafCallback (
     const Eigen::VectorXd& center,
     const Eigen::VectorXd& halflengths,
     int                    position,
     const mesh::Group&     content );

//   /**
//    * @brief Exports the cells of the spacetree into a vtk file.
//    */
//   void doExport ( const Spacetree& spacetree );

private:

   std::string _location;

   std::string _filename;

   size_t _vertexCounter;

   std::list<Eigen::VectorXd> _vertices;

   std::list<std::vector<int> > _cells;

   std::list<int> _cellPositions;

   std::list<int> _cellContents;
};

}} // namespace precice, spacetree

