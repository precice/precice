#include "MeshTest.hpp"
#include "mesh/Vertex.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Quad.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/PropertyContainer.hpp"
#include "mesh/Data.hpp"
#include "utils/Parallel.hpp"
#include "utils/Dimensions.hpp"
#include "tarch/la/WrappedVector.h"
#include "com/MPIDirectCommunication.hpp"
#include "tarch/tests/TestCaseFactory.h"
#include "Eigen/Dense"



registerTest(precice::mesh::tests::MeshTest)

namespace precice {
namespace mesh {
namespace tests {

logging::Logger MeshTest:: _log("precice::mesh::MeshTest");

MeshTest:: MeshTest()
:
  TestCase("mesh::MeshTest")
{}

void MeshTest:: run()
{
  PRECICE_MASTER_ONLY {
    testMethod(testBasicSetup);
    testMethod(testSubIDs);
    testMethod(testComputeState);
    testMethod(testDemonstration);
    testMethod(testBoundingBoxCOG);
  }
# ifndef PRECICE_NO_MPI
  typedef utils::Parallel Par;
  if (Par::getCommunicatorSize() > 3){
    std::vector<int> ranksWanted;
    ranksWanted += 0, 1, 2 , 3;
    MPI_Comm comm = Par::getRestrictedCommunicator(ranksWanted);
    if (Par::getProcessRank() <= 3){
      Par::setGlobalCommunicator(comm);
      testMethod (testDistribution);
      Par::setGlobalCommunicator(Par::getCommunicatorWorld());
    }
  }
# endif // not PRECICE_NO_MPI
}

void MeshTest:: testBasicSetup()
{
  preciceTrace("testBasicSetup()");
  for (int dim=2; dim <= 3; dim++){
    Mesh mesh("MyMesh", dim, false);
    Vertex& v1 = mesh.createVertex(utils::DynVector(dim, 0.0));
    Vertex& v2 = mesh.createVertex(utils::DynVector(dim, 1.0));
    Vertex& v3 = mesh.createVertex(utils::DynVector(dim, 2.0)); // ill-shaped triangle
    Edge& e1 = mesh.createEdge(v1, v2);
    Edge& e2 = mesh.createEdge(v2, v3);
    Edge& e3 = mesh.createEdge(v3, v1);
    if (dim == 3) mesh.createTriangle(e1, e2, e3);
  }
}

void MeshTest:: testSubIDs()
{
  preciceTrace("testSubIDs()");
  for (int dim=2; dim <= 3; dim++){
    Mesh mesh("MyMesh", dim, false);
    Vertex& v0 = mesh.createVertex(utils::DynVector(dim, 0.0));
    Vertex& v1 = mesh.createVertex(utils::DynVector(dim, 1.0));
    PropertyContainer& cont0 = mesh.setSubID("subID0");
    PropertyContainer& cont1 = mesh.setSubID("subID1");
    v0.addParent(cont0);
    v1.addParent(cont0);
    v1.addParent(cont1);
    std::vector<int> properties;
    v0.getProperties(PropertyContainer::INDEX_GEOMETRY_ID, properties);
    validateEquals(properties.size(), 2);
    validateEquals(properties[0], mesh.getID("MyMesh"));
    validateEquals(properties[1], mesh.getID("MyMesh-subID0"));
  }
}

void MeshTest:: testComputeState()
{
  preciceTrace ("testComputeState()");
  bool flipNormals = true;
  { // 2D
    using utils::Vector2D;
    Mesh mesh ( "MyMesh", 2, flipNormals );
    // Create mesh
    Vertex& v1 = mesh.createVertex ( Vector2D(0.0, 0.0) );
    Vertex& v2 = mesh.createVertex ( Vector2D(1.0, 0.0) );
    Vertex& v3 = mesh.createVertex ( Vector2D(1.0, 1.0) );
    //
    //
    // *****
    Edge& e1 = mesh.createEdge ( v1, v2 );
    //     *
    //     *  <---
    // *****
    Edge& e2 = mesh.createEdge ( v2, v3 );
    mesh.computeState();

    // Perform test validations
    using tarch::la::equals;
    validate ( equals(e1.getCenter(), Vector2D(0.5, 0.0)) );
    validate ( equals(e2.getCenter(), Vector2D(1.0, 0.5)) );
    validateNumericalEquals ( e1.getEnclosingRadius(), 0.5 );
    validateNumericalEquals ( e2.getEnclosingRadius(), 0.5 );
    validate ( equals(e1.getNormal(), Vector2D(0.0, 1.0)) );
    validate ( equals(e2.getNormal(), Vector2D(-1.0, 0.0)) );
    validate ( equals(v1.getNormal(), Vector2D(0.0, 1.0)) );
    validate ( equals(v2.getNormal(), Vector2D(-std::sqrt(0.5), std::sqrt(0.5))) );
    validate ( equals(v3.getNormal(), Vector2D(-1.0, 0.0)) );
  }

  { // 3D triangle
    using utils::Vector3D;
    Mesh mesh ( "MyMesh", 3, flipNormals );
    // Create mesh
    Vertex& v1 = mesh.createVertex ( Vector3D(0.0, 0.0, 0.0) );
    Vertex& v2 = mesh.createVertex ( Vector3D(1.0, 0.0, 1.0) );
    Vertex& v3 = mesh.createVertex ( Vector3D(1.0, 1.0, 1.0) );
    Vertex& v4 = mesh.createVertex ( Vector3D(2.0, 0.0, 2.0) );
    Edge& e1 = mesh.createEdge ( v1, v2 );
    Edge& e2 = mesh.createEdge ( v2, v3 );
    Edge& e3 = mesh.createEdge ( v3, v1 );
    Edge& e4 = mesh.createEdge ( v2, v4 );
    Edge& e5 = mesh.createEdge ( v4, v3 );

    //       *
    //     * *
    //   *   *
    // *******
    Triangle& t1 = mesh.createTriangle ( e1, e2, e3 );
    //       *
    //     * * *     <---
    //   *   *   *
    // *************
    Triangle& t2 = mesh.createTriangle ( e4, e5, e2 );
    mesh.computeState();

    // Perform test validations
    using tarch::la::equals;
    validate ( equals(e1.getCenter(), Vector3D(0.5, 0.0, 0.5)) );
    validate ( equals(e2.getCenter(), Vector3D(1.0, 0.5, 1.0)) );
    validate ( equals(e3.getCenter(), Vector3D(0.5, 0.5, 0.5)) );
    validate ( equals(e4.getCenter(), Vector3D(1.5, 0.0, 1.5)) );
    validate ( equals(e5.getCenter(), Vector3D(1.5, 0.5, 1.5)) );
    validateNumericalEquals ( e1.getEnclosingRadius(), std::sqrt(2.0)*0.5 );
    validateNumericalEquals ( e2.getEnclosingRadius(), 0.5 );
    validateNumericalEquals ( e3.getEnclosingRadius(), std::sqrt(3.0)*0.5 );
    validateNumericalEquals ( e4.getEnclosingRadius(), std::sqrt(2.0)*0.5 );
    validateNumericalEquals ( e5.getEnclosingRadius(), std::sqrt(3.0)*0.5 );

    validate ( equals(t1.getCenter(), Vector3D(2.0/3.0, 1.0/3.0, 2.0/3.0)) );
    validate ( equals(t2.getCenter(), Vector3D(4.0/3.0, 1.0/3.0, 4.0/3.0)) );
    validateNumericalEquals ( t1.getEnclosingRadius(), 1.0 );
    validateNumericalEquals ( t2.getEnclosingRadius(), 1.0 );
    Vector3D normal ( 1.0, 0.0, -1.0 );
    normal = normal / tarch::la::norm2(normal);
    validateNumericalEquals ( tarch::la::norm2(normal), 1.0 );
    validate ( equals(t1.getNormal(), normal) );
    validate ( equals(t2.getNormal(), normal) );

    validate ( equals(e1.getNormal(), normal) );
    validate ( equals(e2.getNormal(), normal) );
    validate ( equals(e3.getNormal(), normal) );
    validate ( equals(e4.getNormal(), normal) );
    validate ( equals(e5.getNormal(), normal) );
    validate ( equals(v1.getNormal(), normal) );
    validate ( equals(v2.getNormal(), normal) );
    validate ( equals(v3.getNormal(), normal) );
    validate ( equals(v4.getNormal(), normal) );
  }

  { // 3D quad
    using utils::Vector3D;
    Mesh mesh("MyMesh", 3, flipNormals);
    // Create mesh (Two rectangles with a common edge at z-axis. One extends in
    // x-y-plane (2 long), the other in y-z-plane (1 long))
    Vertex& v0 = mesh.createVertex(Vector3D(0.0, 0.0, 0.0));
    Vertex& v1 = mesh.createVertex(Vector3D(2.0, 0.0, 0.0));
    Vertex& v2 = mesh.createVertex(Vector3D(2.0, 1.0, 0.0));
    Vertex& v3 = mesh.createVertex(Vector3D(0.0, 1.0, 0.0));
    Vertex& v4 = mesh.createVertex(Vector3D(0.0, 0.0, 1.0));
    Vertex& v5 = mesh.createVertex(Vector3D(0.0, 1.0, 1.0));
    Edge& e0 = mesh.createEdge(v0, v1);
    Edge& e1 = mesh.createEdge(v1, v2);
    Edge& e2 = mesh.createEdge(v2, v3);
    Edge& e3 = mesh.createEdge(v3, v0);
    Edge& e4 = mesh.createEdge(v0, v4);
    Edge& e5 = mesh.createEdge(v4, v5);
    Edge& e6 = mesh.createEdge(v5, v3);

    Quad& q0 = mesh.createQuad(e0, e1, e2, e3); // in x-y-plane
    Quad& q1 = mesh.createQuad(e4, e5, e6, e3); // in z-y-plane
    mesh.computeState();

    // Perform test validations
    using tarch::la::equals;

    validate(equals(e0.getCenter(), Vector3D(1.0, 0.0, 0.0)));
    validate(equals(e1.getCenter(), Vector3D(2.0, 0.5, 0.0)));
    validate(equals(e2.getCenter(), Vector3D(1.0, 1.0, 0.0)));
    validate(equals(e3.getCenter(), Vector3D(0.0, 0.5, 0.0)));
    validate(equals(e4.getCenter(), Vector3D(0.0, 0.0, 0.5)));
    validate(equals(e5.getCenter(), Vector3D(0.0, 0.5, 1.0)));
    validate(equals(e6.getCenter(), Vector3D(0.0, 1.0, 0.5)));

    validateNumericalEquals(e0.getEnclosingRadius(), 1.0);
    validateNumericalEquals(e1.getEnclosingRadius(), 0.5);
    validateNumericalEquals(e2.getEnclosingRadius(), 1.0);
    validateNumericalEquals(e3.getEnclosingRadius(), 0.5);
    validateNumericalEquals(e4.getEnclosingRadius(), 0.5);
    validateNumericalEquals(e5.getEnclosingRadius(), 0.5);
    validateNumericalEquals(e6.getEnclosingRadius(), 0.5);

    validate(equals(q0.getCenter(), Vector3D(1.0, 0.5, 0.0)));
    validate(equals(q1.getCenter(), Vector3D(0.0, 0.5, 0.5)));

    validateNumericalEquals(q0.getEnclosingRadius(), sqrt(1.25));
    validateNumericalEquals(q1.getEnclosingRadius(), sqrt(0.5));

    Vector3D normal0(0.0, 0.0, -1.0);
    validateWithParams1(equals(q0.getNormal(), normal0), q0.getNormal());
    Vector3D normal1(1.0, 0.0, 0.0);
    validateWithParams1(equals(q1.getNormal(), normal1), q1.getNormal());

    validate(equals(e0.getNormal(), normal0));
    validate(equals(e1.getNormal(), normal0));
    validate(equals(e2.getNormal(), normal0));
    Vector3D normal0and1(2.0*normal0 + normal1);
    normal0and1 /= norm2(normal0and1);
    validateWithParams2(equals(e3.getNormal(), normal0and1), e3.getNormal(), normal0and1);
    validate(equals(e4.getNormal(), normal1));
    validate(equals(e5.getNormal(), normal1));
    validate(equals(e6.getNormal(), normal1));

    validate(equals(v0.getNormal(), normal0and1));
    validate(equals(v1.getNormal(), normal0));
    validate(equals(v2.getNormal(), normal0));
    validate(equals(v3.getNormal(), normal0and1));
    validate(equals(v4.getNormal(), normal1));
    validate(equals(v5.getNormal(), normal1));
  }
}

void MeshTest::testBoundingBoxCOG()
{
  {
    utils::DynVector coords0(2);
    utils::DynVector coords1(2);
    utils::DynVector coords2(2);

    coords0 =  2.0, 0.0;
    coords1 = -1.0, 4.0;
    coords2 =  0.0, 1.0;

    Mesh mesh ("2D Testmesh", 2, false );
    mesh.createVertex(coords0);
    mesh.createVertex(coords1);
    mesh.createVertex(coords2);

    mesh.computeState();

    Mesh::BoundingBox bBox = mesh.getBoundingBox();
    auto cog = mesh.getCOG();

    Mesh::BoundingBox referenceBox =  { {-1.0, 2.0},
                                        { 0.0, 4.0} };

    std::vector<double> referenceCOG =  { 0.5, 2.0 };

    validateEquals(bBox.size(), 2);
    validateEquals(cog.size(), 2);

    for (size_t d = 0; d < bBox.size(); d++) {
      validateNumericalEquals(referenceBox[d].first, bBox[d].first);
      validateNumericalEquals(referenceBox[d].second, bBox[d].second);
    }
    for (size_t d = 0; d < cog.size(); d++) {
      validateNumericalEquals(referenceCOG[d], cog[d]);
    }
  }

  {
    utils::DynVector coords0(3);
    utils::DynVector coords1(3);
    utils::DynVector coords2(3);
    utils::DynVector coords3(3);

    coords0 =  2.0, 0.0, -3.0;
    coords1 = -1.0, 4.0,  8.0;
    coords2 =  0.0, 1.0, -2.0;
    coords3 =  3.5, 2.0, -2.0;

    Mesh mesh ("3D Testmesh", 3, false );
    mesh.createVertex(coords0);
    mesh.createVertex(coords1);
    mesh.createVertex(coords2);
    mesh.createVertex(coords3);

    mesh.computeState();

    Mesh::BoundingBox bBox = mesh.getBoundingBox();
    auto cog = mesh.getCOG();

    Mesh::BoundingBox referenceBox =  { {-1.0, 3.5},
                                        { 0.0, 4.0},
                                        {-3.0, 8.0} };

    std::vector<double> referenceCOG =  { 1.25, 2.0, 2.5 };

    validateEquals(bBox.size(), 3);
    validateEquals(cog.size(), 3);

    for (size_t d = 0; d < bBox.size(); d++) {
      validateNumericalEquals(referenceBox[d].first, bBox[d].first);
      validateNumericalEquals(referenceBox[d].second, bBox[d].second);
    }
    for (size_t d = 0; d < cog.size(); d++) {
      validateNumericalEquals(referenceCOG[d], cog[d]);
    }
  }
}


void MeshTest:: testDemonstration ()
{
  preciceTrace ( "testDemonstration()" );
  for ( int dim=2; dim <= 3; dim++ ){
    preciceDebug ( "dim = " << dim );
    // Create mesh object
    std::string meshName ( "MyMesh" );
    bool flipNormals = false; // The normals of triangles, edges, vertices
    Mesh mesh ( meshName, dim, flipNormals );
    int geometryID = mesh.getProperty<int> ( mesh.INDEX_GEOMETRY_ID );

    // Validate mesh object state
    validateEquals(mesh.getName(), meshName);
    validateEquals(mesh.isFlipNormals(), flipNormals);

    // Create mesh vertices
    utils::DynVector coords0(dim);
    utils::DynVector coords1(dim);
    utils::DynVector coords2(dim);
    if (dim == 2){
      coords0 = 0.0, 0.0;
      coords1 = 1.0, 0.0;
      coords2 = 0.0, 1.0;
    }
    else {
      coords0 = 0.0, 0.0, 0.0;
      coords1 = 1.0, 0.0, 0.0;
      coords2 = 0.0, 0.0, 1.0;
    }
    Vertex& v0 = mesh.createVertex(coords0);
    Vertex& v1 = mesh.createVertex(coords1);
    Vertex& v2 = mesh.createVertex(coords2);

    // Validate mesh vertices state
    // This is the preferred way to iterate over elements in a mesh, it hides
    // the details of the vertex container in class Mesh.
    size_t index = 0;
    for (Vertex& vertex : mesh.vertices()) {
      if ( index == 0 ) {
        validateEquals ( vertex.getID(), v0.getID() );
      }
      else if ( index == 1 ) {
        validateEquals ( vertex.getID(), v1.getID() );
      }
      else if ( index == 2 ) {
        validateEquals ( vertex.getID(), v2.getID() );
      }
      else {
        validate ( false );
      }
      index++;
    }

    // Create mesh edges
    Edge& e0 = mesh.createEdge ( v0, v1 );
    Edge& e1 = mesh.createEdge ( v1, v2 );
    Edge& e2 = mesh.createEdge ( v2, v0 );

    // Validate mesh edges state
    index = 0;
    for (Edge & edge : mesh.edges()) {
      if ( index == 0 ) {
        validateEquals ( edge.getID(), e0.getID() );
      }
      else if ( index == 1 ) {
        validateEquals ( edge.getID(), e1.getID() );
      }
      else if ( index == 2 ) {
        validateEquals ( edge.getID(), e2.getID() );
      }
      else {
        validate ( false );
      }
      index++;
    }

    Triangle* t = nullptr;
    if ( dim == 3 ){
      // Create triangle
      t = & mesh.createTriangle ( e0, e1, e2 );

      // Validate mesh triangle
      validateEquals ( (*mesh.triangles().begin()).getID(), t->getID() );
    }

    // Create vertex data
    std::string dataName ( "MyData" );
    int dataDimensions = dim;
    // Add a data set to the mesh. Every data value is associated to a vertex in
    // the mesh via the vertex ID. The data values are created when
    // Mesh::computeState() is called by the mesh holding the data.
    PtrData data = mesh.createData ( dataName, dataDimensions );

    // Validate data state
    validateEquals ( data->getName(), dataName );
    validateEquals ( data->getDimensions(), dataDimensions );


    // Validate state of mesh with data
    validateEquals ( mesh.data().size(), 1 );
    validateEquals ( mesh.data()[0]->getName(), dataName );

    // Create sub-id
    std::string nameSubIDPrefix ( "sub-id" );
    PropertyContainer & cont = mesh.setSubID ( nameSubIDPrefix );
    int subID = cont.getFreePropertyID();
    cont.setProperty ( cont.INDEX_GEOMETRY_ID, subID );

    // Add sub-id to selected mesh elements
    v0.addParent ( cont );
    v1.addParent ( cont );
    e0.addParent ( cont );

    // Validate geometry IDs
    std::vector<int> geometryIDs;
    v0.getProperties ( v0.INDEX_GEOMETRY_ID, geometryIDs );
    validateEquals ( geometryIDs.size(), 2 );
    validate ( utils::contained(geometryID, geometryIDs) );
    validate ( utils::contained(subID, geometryIDs) );
    geometryIDs.clear();
    v1.getProperties ( v1.INDEX_GEOMETRY_ID, geometryIDs );
    validateEquals ( geometryIDs.size(), 2 );
    validate ( utils::contained(geometryID, geometryIDs) );
    validate ( utils::contained(subID, geometryIDs) );
    geometryIDs.clear();
    v2.getProperties ( v1.INDEX_GEOMETRY_ID, geometryIDs );
    validateEquals ( geometryIDs.size(), 1 );
    validate ( utils::contained(geometryID, geometryIDs) );
    geometryIDs.clear();
    e0.getProperties ( e0.INDEX_GEOMETRY_ID, geometryIDs );
    validateEquals ( geometryIDs.size(), 2 );
    validate ( utils::contained(geometryID, geometryIDs) );
    validate ( utils::contained(subID, geometryIDs) );
    geometryIDs.clear();
    e1.getProperties ( e1.INDEX_GEOMETRY_ID, geometryIDs );
    validateEquals ( geometryIDs.size(), 1 );
    validate ( utils::contained(geometryID, geometryIDs) );
    geometryIDs.clear();
    e2.getProperties ( e2.INDEX_GEOMETRY_ID, geometryIDs );
    validateEquals ( geometryIDs.size(), 1 );
    validate ( utils::contained(geometryID, geometryIDs) );
    if ( dim == 3 ){
      geometryIDs.clear();
      t->getProperties ( t->INDEX_GEOMETRY_ID, geometryIDs );
      validateEquals ( geometryIDs.size(), 1 );
      validate ( utils::contained(geometryID, geometryIDs) );
    }
    validateEquals ( mesh.getNameIDPairs().size(), 2 );
    validateEquals ( mesh.getNameIDPairs().count("MyMesh"), 1 );
    validateEquals ( mesh.getNameIDPairs().count("MyMesh-sub-id"), 1 );

    // Compute the state of the mesh elements (vertices, edges, triangles)
    mesh.computeState();

    // Allocate memory for the data values of set data. Before data value access
    // leads to assertions.
    mesh.allocateDataValues();

    // Access data values
    Eigen::VectorXd& dataValues = data->values();
    validateEquals ( dataValues.size(), 3 * dim );
    validateEquals ( v0.getID(), 0 );
    using tarch::la::slice;
    Eigen::VectorXd value = Eigen::VectorXd::Zero(dim);
    for ( int i=0; i < dim; i++ ){
      value(i) = dataValues(v0.getID() * dim + i);
    }
  }
}

void MeshTest:: testDistribution()
{
  preciceTrace ("testDistribution()");
# ifndef PRECICE_NO_MPI
  assertion ( utils::Parallel::getCommunicatorSize() == 4 );

  com::Communication::SharedPointer masterSlaveCom =
      com::Communication::SharedPointer(new com::MPIDirectCommunication());
  utils::MasterSlave::_communication = masterSlaveCom;

  if (utils::Parallel::getProcessRank() == 0){ //Master
    utils::Parallel::splitCommunicator( "Master" );
    utils::MasterSlave::_rank = 0;
    utils::MasterSlave::_size = 4;
    utils::MasterSlave::_slaveMode = false;
    utils::MasterSlave::_masterMode = true;
    masterSlaveCom->acceptConnection("Master", "Slaves", 0, 1);
    masterSlaveCom->setRankOffset(1);
  }
  else if(utils::Parallel::getProcessRank() == 1){
    utils::Parallel::splitCommunicator( "Slaves" );
    utils::MasterSlave::_rank = 1;
    utils::MasterSlave::_size = 4;
    utils::MasterSlave::_slaveMode = true;
    utils::MasterSlave::_masterMode = false;
    masterSlaveCom->requestConnection("Master", "Slaves", 0, 3);
  }
  else if(utils::Parallel::getProcessRank() == 2){
    utils::Parallel::splitCommunicator( "Slaves");
    utils::MasterSlave::_rank = 2;
    utils::MasterSlave::_size = 4;
    utils::MasterSlave::_slaveMode = true;
    utils::MasterSlave::_masterMode = false;
    masterSlaveCom->requestConnection("Master", "Slaves", 1, 3);
  }
  else if(utils::Parallel::getProcessRank() == 3){
    utils::Parallel::splitCommunicator( "Slaves");
    utils::MasterSlave::_rank = 3;
    utils::MasterSlave::_size = 4;
    utils::MasterSlave::_slaveMode = true;
    utils::MasterSlave::_masterMode = false;
    masterSlaveCom->requestConnection("Master", "Slaves", 2, 3);
  }

  // Create mesh object
  std::string meshName ( "MyMesh" );
  int dim = 2;
  bool flipNormals = false; // The normals of triangles, edges, vertices
  Mesh mesh ( meshName, dim, flipNormals );

  if (utils::Parallel::getProcessRank() == 0) { //Master
    utils::DynVector position(dim);
    assignList(position) = 0.0, 0.0;
    mesh.createVertex(position);
    assignList(position) = 1.0, 0.0;
    mesh.createVertex(position);
  } else if (utils::Parallel::getProcessRank() == 1) { //Slave1
    utils::DynVector position(dim);
    assignList(position) = 1.0, 0.0;
    mesh.createVertex(position);
  } else if (utils::Parallel::getProcessRank() == 2) { //Slave2
  } else if (utils::Parallel::getProcessRank() == 3) { //Slave3
    utils::DynVector position(dim);
    assignList(position) = 1.0, 0.0;
    mesh.createVertex(position);
    assignList(position) = 2.0, 0.0;
    mesh.createVertex(position);
  }

  if (utils::Parallel::getProcessRank() == 0) { //Master
    mesh.getVertexDistribution()[0].push_back(0);
    mesh.getVertexDistribution()[0].push_back(1);
    mesh.getVertexDistribution()[1].push_back(1);
    mesh.getVertexDistribution()[3].push_back(1);
    mesh.getVertexDistribution()[3].push_back(2);
    mesh.setGlobalNumberOfVertices(3);
  }

  mesh.computeDistribution();

  if (utils::Parallel::getProcessRank() == 0) { //Master
    validate(mesh.getGlobalNumberOfVertices()==3);
    validate(mesh.getVertexOffsets().size()==4);
    validate(mesh.getVertexOffsets()[0]==2);
    validate(mesh.getVertexOffsets()[1]==3);
    validate(mesh.getVertexOffsets()[2]==3);
    validate(mesh.getVertexOffsets()[3]==5);
    validate(mesh.vertices()[0].getGlobalIndex()==0);
    validate(mesh.vertices()[1].getGlobalIndex()==1);
    validate(mesh.vertices()[0].isOwner()==true);
    validate(mesh.vertices()[1].isOwner()==true);
  } else if (utils::Parallel::getProcessRank() == 1) { //Slave1
    validate(mesh.getGlobalNumberOfVertices()==3);
    validate(mesh.getVertexOffsets().size()==4);
    validate(mesh.getVertexOffsets()[0]==2);
    validate(mesh.getVertexOffsets()[1]==3);
    validate(mesh.getVertexOffsets()[2]==3);
    validate(mesh.getVertexOffsets()[3]==5);
    validate(mesh.vertices()[0].getGlobalIndex()==1);
    validate(mesh.vertices()[0].isOwner()==false);
  } else if (utils::Parallel::getProcessRank() == 2) { //Slave2
    validate(mesh.getGlobalNumberOfVertices()==3);
    validate(mesh.getVertexOffsets().size()==4);
    validate(mesh.getVertexOffsets()[0]==2);
    validate(mesh.getVertexOffsets()[1]==3);
    validate(mesh.getVertexOffsets()[2]==3);
    validate(mesh.getVertexOffsets()[3]==5);
  } else if (utils::Parallel::getProcessRank() == 3) { //Slave3
    validate(mesh.getGlobalNumberOfVertices()==3);
    validate(mesh.getVertexOffsets().size()==4);
    validate(mesh.getVertexOffsets()[0]==2);
    validate(mesh.getVertexOffsets()[1]==3);
    validate(mesh.getVertexOffsets()[2]==3);
    validate(mesh.getVertexOffsets()[3]==5);
    validate(mesh.vertices()[0].getGlobalIndex()==1);
    validate(mesh.vertices()[1].getGlobalIndex()==2);
    validate(mesh.vertices()[0].isOwner()==false);
    validate(mesh.vertices()[1].isOwner()==true);
  }

  utils::MasterSlave::_slaveMode = false;
  utils::MasterSlave::_masterMode = false;
  utils::Parallel::synchronizeProcesses();
  utils::Parallel::clearGroups();
  utils::MasterSlave::_communication = nullptr;
#endif // PRECICE_NO_MPI
}

}}} // namespace precice, mesh, tests
