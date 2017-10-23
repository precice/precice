#ifndef PRECICE_QUERY_FINDCLOSESTTEST_HPP_
#define PRECICE_QUERY_FINDCLOSESTTEST_HPP_

#include "tarch/tests/TestCase.h"
#include "logging/Logger.hpp"

namespace precice {
namespace query {
namespace tests {

class FindClosestTest : public tarch::tests::TestCase
{
public:

  /**
   * @brief Constructor.
   */
	FindClosestTest();

	/**
	 * @brief Prepares tests. Empty.
	 */
	virtual void setUp() {}

	/**
	 * @brief Runs all tests.
	 */
  virtual void run();

private:

    /**
     * @brief Tests finding shortest distances with Vertex objects only.
     */
    void testFindClosestDistanceToVertices();

    /**
     * @brief Tests sign of shortest distances with Vertex objects only.
     */
    void testSignOfShortestDistance();

    /**
     * @brief Tests if abs value of shortest distance is independent from sign.
     */
    void testIndependenceOfSignOfShortestDistance();

    /**
     * @brief Tests finding shortests distances with Vertex and Edge objects.
     */
    void testFindClosestDistanceToEdges();

    void testFindClosestDistanceToEdges3D();

    /**
     * @brief Tests finding shortests distances with Vertex, Edge and Triangle object.
     */
    void testFindClosestDistanceToTriangles();

    void testFindClosestDistanceToTrianglesAndVertices();

    /**
     * @brief Tests finding the correct shortest distance to quadrangles.
     */
    //void testFindClosestDistanceToQuads();

    /**
     * @brief Tests and visualizes queries on meshes with different meshIDs
     *
     * The mesh consists of one triangle with 2/3 vertices.
     * The vertex coordinates are created in a loop, where each vertex coordinate
     * is set either 1.0, if the dimension of the coordinate equals the vertex's
     * (local) index, or 0.0 otherwise.
     */
    void testMultipleMeshIDs();

    void testWeigthsOfVertices();

private:

   // @brief Logging device.
   static logging::Logger _log;
};

}}} // namespace precice, query, tests

#endif /* PRECICE_QUERY_FINDCLOSESTTEST_HPP_ */
