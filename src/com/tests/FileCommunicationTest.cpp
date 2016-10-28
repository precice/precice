#include "FileCommunicationTest.hpp"
#include "com/FileCommunication.hpp"
#include "utils/Globals.hpp"
#include "utils/Parallel.hpp"
#include "math/math.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::com::tests::FileCommunicationTest)

namespace precice {
namespace com {
namespace tests {

logging::Logger FileCommunicationTest::
  _log ( "precice::com::tests::FileCommunicationTest" );

FileCommunicationTest:: FileCommunicationTest()
:
  TestCase ( "precice::com::tests::FileCommunicationTest" )
{}

void FileCommunicationTest:: run()
{
  preciceTrace ( "run()" );
  typedef utils::Parallel Par;
  if ( Par::getCommunicatorSize() >= 2 ){
    if ( Par::getProcessRank() < 2 ){
      //TODO: are not working on Benjamin's laptop
      //testMethod ( testSimpleSendReceive );
      //testMethod ( testMultipleExchanges );
    }
  }
}

void FileCommunicationTest:: testSimpleSendReceive()
{
  preciceTrace ( "testSimpleSendReceive()" );
  int rank = utils::Parallel::getProcessRank();
  bool binaryMode = false;
  FileCommunication comTxt ( binaryMode, "" );
  validateEquals ( comTxt.isConnected(), false );
  std::string requester ( "FileCommunicationTest-testSimpleSendReceive-Requester-txt" );
  std::string acceptor ( "FileCommunicationTest-testSimpleSendReceive-Acceptor-txt" );
  Eigen::Vector3d doubleVector = Eigen::Vector3d::Constant(0.1234567890123456);
  Eigen::Matrix<int, 5, 1> intVector;
  intVector << 1, 2, 3, 4, 5;
  // int intVector[] = { 1, 2, 3, 4, 5 };
  int integer = 1;
  double dbl = 1.0;
  bool boolean = true;
  if ( rank == 0 ){
    comTxt.requestConnection ( requester, acceptor, 0, 1 );
    validate ( comTxt.isConnected() );
    comTxt.startSendPackage ( 0 );
    comTxt.send ( std::string("rank0"), 0 );
    comTxt.send ( doubleVector.data(), doubleVector.size(), 0 );
    comTxt.send ( intVector.data(), 5, 0 );
    comTxt.send ( dbl, 0 );
    comTxt.send ( integer, 0 );
    comTxt.send ( boolean, 0 );
    comTxt.finishSendPackage ();
  }
  else {
    assertion ( rank == 1, rank );
    comTxt.acceptConnection ( requester, acceptor, 0, 1 );
    validate ( comTxt.isConnected() );
    comTxt.startReceivePackage ( 0 );
    std::string message ( "" ) ;
    comTxt.receive ( message, 0 );
    double rawVec[3];
    // assign(wrap<3>(rawVec)) = 0.0;
    Eigen::Map<Eigen::Vector3d>(rawVec).setConstant(0);
    comTxt.receive ( rawVec, 3, 0 );
    int intArrayReceived[] = { 0, 0, 0, 0, 0 };
    comTxt.receive ( intArrayReceived, 5, 0 );
    double dblReceived = 0.0;
    comTxt.receive ( dblReceived, 0 );
    int integerReceived = 0;
    comTxt.receive ( integerReceived, 0 );
    bool booleanReceived = false;
    comTxt.receive ( booleanReceived, 0 );
    comTxt.finishReceivePackage ();
    validateEquals ( message, std::string("rank0") );
    validate ( math::equals(Eigen::Map<Eigen::Vector3d>(rawVec), doubleVector));
    validate ( math::equals(Eigen::Map<Eigen::Matrix<int, 5, 1>>(intArrayReceived), intVector));
    validate ( math::equals(dblReceived, dbl) );
    validateEquals ( integerReceived, integer );
    validateEquals ( booleanReceived, boolean );
  }
  comTxt.closeConnection();
  validate ( not comTxt.isConnected() );

  binaryMode = true;
  FileCommunication comBin ( binaryMode, "" );
  validateEquals ( comBin.isConnected(), false );
  requester = "FileCommunicationTest-testSimpleSendReceive-Requester-bin";
  acceptor = "FileCommunicationTest-testSimpleSendReceive-Acceptor-bin";
  if ( rank == 0 ){
    comBin.requestConnection ( requester, acceptor, 0, 1 );
    validate ( comBin.isConnected() );
    comBin.startSendPackage ( 0 );
    comBin.send ( std::string("rank0"), 0 );
    comBin.send ( doubleVector.data(), 3, 0 );
    comBin.send ( intVector.data(), 5, 0 );
    comBin.send ( dbl, 0 );
    comBin.send ( integer, 0 );
    comBin.send ( boolean, 0 );
    comBin.finishSendPackage();
  }
  else {
    assertion ( rank == 1, rank );
    comBin.acceptConnection ( requester, acceptor, 0, 1 );
    validate ( comBin.isConnected() );
    comBin.startReceivePackage ( 0 );
    std::string message ( "" ) ;
    comBin.receive ( message, 0 );
    double rawVec[3];
    Eigen::Map<Eigen::Vector3d>(rawVec).setConstant(0);
    comBin.receive ( rawVec, 3, 0 );
    int intArrayReceived[] = { 0, 0, 0, 0, 0 };
    comBin.receive ( intArrayReceived, 5, 0 );
    double dblReceived = 0.0;
    comBin.receive ( dblReceived, 0 );
    int integerReceived = 0;
    comBin.receive ( integerReceived, 0 );
    bool booleanReceived = false;
    comBin.receive ( booleanReceived, 0 );
    comBin.finishReceivePackage ();
    validateEquals ( message, std::string("rank0") );
    validate ( math::equals(Eigen::Map<Eigen::Vector3d>(rawVec), doubleVector));
    validate ( math::equals(Eigen::Map<Eigen::Matrix<int, 5, 1>>(intArrayReceived), intVector));
    validate ( math::equals(dblReceived, dbl) );
    validateEquals ( integerReceived, integer );
    validateEquals ( booleanReceived, boolean );
  }
  comBin.closeConnection ();
  validate ( not comBin.isConnected() );
}

void FileCommunicationTest:: testMultipleExchanges()
{
  preciceTrace ( "testMultipleExchanges()" );
  int rank = utils::Parallel::getProcessRank();
  bool binaryMode = false;
  FileCommunication com ( binaryMode, "" );
  validate ( not com.isConnected() );
  std::string requester ( "FileCommunicationTest-testMultipleExchanges-Requester-txt" );
  std::string acceptor ( "FileCommunicationTest-testMultipleExchanges-Acceptor-txt" );
  int value0 = 1;
  int value1 = 2;
  int value2 = 3;
  if ( rank == 0 ){
    com.requestConnection ( acceptor, requester, 0, 1 );
    validate ( com.isConnected() );

    com.startSendPackage ( 0 );
    com.send ( value0, 0 );
    com.finishSendPackage ();

    com.startReceivePackage ( 0 );
    int number = 0;
    com.receive ( number, 0 );
    com.finishReceivePackage ();
    validateEquals ( number, value1 );

    com.startSendPackage ( 0 );
    com.send ( value2, 0 );
    com.finishSendPackage ();

    com.startSendPackage ( 0 );
    com.send ( value0, 0 );
    com.finishSendPackage ();

    com.startSendPackage ( 0 );
    com.send ( value1, 0 );
    com.finishSendPackage ();

    com.startReceivePackage ( 0 );
    com.receive ( number, 0 );
    com.finishReceivePackage ();
    validateEquals ( number, value2 );

    com.closeConnection ();
    validate ( not com.isConnected() );
  }
  else {
    assertion ( rank == 1, rank );
    com.acceptConnection ( acceptor, requester, 0, 1 );
    validate ( com.isConnected() );

    com.startReceivePackage ( 0 );
    int number = 0;
    com.receive ( number, 0 );
    com.finishReceivePackage ();
    validateEquals ( number, value0 );

    com.startSendPackage ( 0 );
    com.send ( value1, 0 );
    com.finishSendPackage ();

    com.startReceivePackage ( 0 );
    com.receive ( number, 0 );
    com.finishReceivePackage ();
    validateEquals ( number, value2 );

    com.startReceivePackage ( 0 );
    com.receive ( number, 0 );
    com.finishReceivePackage ();
    validateEquals ( number, value0 );

    com.startReceivePackage ( 0 );
    com.receive ( number, 0 );
    com.finishReceivePackage ();
    validateEquals ( number, value1 );

    com.startSendPackage ( 0 );
    com.send ( value2, 0 );
    com.finishSendPackage ();

    com.closeConnection ();
    validate ( not com.isConnected() );
  }
}

}}} // namespace precice, com, tests
