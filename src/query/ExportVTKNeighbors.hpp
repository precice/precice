#ifndef PRECICE_QUERY_EXPORTVTKNEIGHBORS_HPP_
#define PRECICE_QUERY_EXPORTVTKNEIGHBORS_HPP_

#include "FindClosest.hpp"
#include <string>

namespace precice {
namespace query {

/*
 * @brief Visualizes results of FindClosestOnMeshVisitor searches
 *
 * The visualization consists of lines connecting search points with
 * closest points on found visitables. The lines have assigned a value
 * distance, that has the value of the distance to the found visitables.
 *
 * @author Bernhard Gatzhammer
 */
class ExportVTKNeighbors
{
public:

   /**
    * @brief Adds a neighborhood relation already computed.
    *
    * @param[in] queryPoint Point whose next neighbor was determined.
    * @param[in] neighborPoint ClosestNeighbor of queryPoint.
    */
   void addNeighbors (
     const Eigen::VectorXd&  queryPoint,
     const ClosestElement&   closestNeighbor );

   /**
    * @brief Visualizes results of performed searches into vtk file
    */
   void exportNeighbors ( const std::string& filename );

   /**
    * @brief Resets results of performed searches
    */
   void resetElements() { _neighbors.clear(); };

private:

   /// Results of performed searches
  std::vector<std::pair<Eigen::VectorXd, query::ClosestElement>> _neighbors;
};

}} // namespace precice, query

#endif /* PRECICE_QUERY_EXPORTVTKNEIGHBORS_HPP_ */
