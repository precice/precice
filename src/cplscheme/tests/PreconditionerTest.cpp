
#include "PreconditionerTest.h"
#include "cplscheme/impl/ResidualPreconditioner.hpp"
#include "cplscheme/impl/ResidualSumPreconditioner.hpp"
#include "cplscheme/impl/ValuePreconditioner.hpp"
#include "cplscheme/impl/ConstantPreconditioner.hpp"
#include "cplscheme/impl/SharedPointer.hpp"
#include <Eigen/Dense>
#include "tarch/la/DynamicMatrix.h"
#include "utils/MasterSlave.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "utils/Parallel.hpp"

#include "tarch/tests/TestCaseFactory.h"

registerTest(precice::cplscheme::tests::PreconditionerTest)

namespace precice {
namespace cplscheme {
namespace tests {

logging::Logger PreconditionerTest::
   _log ( "precice::cplscheme::tests::PreconditionerTest" );

PreconditionerTest::PreconditionerTest ()
:
  TestCase ("precice::cplscheme::PreconditionerTest"),
  _data(),
  _res(),
  _compareDataRes(),
  _compareDataResSum(),
  _compareDataValue(),
  _compareDataConstant()
{}

void PreconditionerTest::run ()
{
  preciceTrace ( "run" );

#ifndef PRECICE_NO_MPI
  typedef utils::Parallel Par;
  if (Par::getCommunicatorSize() > 3){
    std::vector<int> ranksWanted;
    ranksWanted += 0, 1, 2 , 3;
    MPI_Comm comm = Par::getRestrictedCommunicator(ranksWanted);
    if (Par::getProcessRank() <= 3){
      Par::setGlobalCommunicator(comm);
      testMethod ( testParallelMatrixScaling );
      Par::setGlobalCommunicator(Par::getCommunicatorWorld());
    }
  }
#endif

  PRECICE_MASTER_ONLY {
    testMethod (testResPreconditioner);
    testMethod (testResSumPreconditioner);
    testMethod (testValuePreconditioner);
    testMethod (testConstPreconditioner);
  }
}


void PreconditionerTest::setUp()
{
  _data.append(1.0);
  _data.append(2.0);
  _data.append(3.0);
  _data.append(4.0);
  _data.append(5.0);
  _data.append(6.0);
  _data.append(7.0);
  _data.append(8.0);
  _res.append(0.1);
  _res.append(0.1);
  _res.append(0.001);
  _res.append(0.001);
  _res.append(0.001);
  _res.append(0.001);
  _res.append(10.0);
  _res.append(20.0);
  _compareDataRes.append(7.07106781186547372897e+00);
  _compareDataRes.append(1.41421356237309474579e+01);
  _compareDataRes.append(1.50000000000000000000e+03);
  _compareDataRes.append(2.00000000000000000000e+03);
  _compareDataRes.append(2.50000000000000000000e+03);
  _compareDataRes.append(3.00000000000000000000e+03);
  _compareDataRes.append(3.13049516849970566046e-01);
  _compareDataRes.append(3.57770876399966353265e-01);
  _compareDataResSum.append(7.90585229434499154877e+01);
  _compareDataResSum.append(1.58117045886899830975e+02);
  _compareDataResSum.append(1.67708453051717078779e+04);
  _compareDataResSum.append(2.23611270735622783832e+04);
  _compareDataResSum.append(2.79514088419528488885e+04);
  _compareDataResSum.append(3.35416906103434157558e+04);
  _compareDataResSum.append(3.50007001329973377324e+00);
  _compareDataResSum.append(4.00008001519969536020e+00);
  _compareDataValue.append(4.47213595499957927704e-01);
  _compareDataValue.append(8.94427190999915855407e-01);
  _compareDataValue.append(3.23498319610315276940e-01);
  _compareDataValue.append(4.31331092813753647075e-01);
  _compareDataValue.append(5.39163866017192239255e-01);
  _compareDataValue.append(6.46996639220630553879e-01);
  _compareDataValue.append(6.58504607868518165859e-01);
  _compareDataValue.append(7.52576694706877713514e-01);
  _compareDataConstant.append(1.00000000000000002082e-03);
  _compareDataConstant.append(2.00000000000000004163e-03);
  _compareDataConstant.append(1.49999999999999977796e+00);
  _compareDataConstant.append(1.99999999999999955591e+00);
  _compareDataConstant.append(2.50000000000000044409e+00);
  _compareDataConstant.append(2.99999999999999955591e+00);
  _compareDataConstant.append(6.99999999999999883585e+05);
  _compareDataConstant.append(7.99999999999999650754e+05);

}

void PreconditionerTest::testResPreconditioner ()
{
  preciceTrace("testResPreconditioner()");
  std::vector<int> dims;
  dims.push_back(1);
  dims.push_back(2);
  dims.push_back(1);

  impl::ResidualPreconditioner precond(dims);

  int numberOfRows = 8;
  precond.initialize(numberOfRows);
  DataValues backup = _data;

  preciceDebug("New iteration");
  //should change
  precond.update(false, _data, _res);
  validate(precond.requireNewQR());
  precond.newQRfulfilled();
  precond.apply(_data);
  validateVector(_data, _compareDataRes);
  precond.revert(_data);
  validateVector(_data, backup);

  preciceDebug("New timestep");
  //should not change weights
  precond.update(true, _data, _res*10);
  validate(not precond.requireNewQR());
  precond.apply(_data);
  validateVector(_data, _compareDataRes);
  precond.revert(_data);
  validateVector(_data, backup);
}

void PreconditionerTest::testResSumPreconditioner ()
{
  preciceTrace("testResSumPreconditioner()");
  std::vector<int> dims;
  dims.push_back(1);
  dims.push_back(2);
  dims.push_back(1);

  impl::ResidualSumPreconditioner precond(dims);

  int numberOfRows = 8;
  precond.initialize(numberOfRows);
  DataValues backup = _data;

  preciceDebug("New iteration");
  //should change, update twice to really test the summation
  precond.update(false, _data, _res);
  precond.update(false, _data, _res*2);
  validate(precond.requireNewQR());
  precond.newQRfulfilled();
  precond.apply(_data);
  validateVector(_data, _compareDataResSum);

  precond.revert(_data);
  validateVector(_data, backup);

  preciceDebug("New timestep");
  //should not change weights
  precond.update(true, _data, _res*10);
  validate(not precond.requireNewQR());
  precond.apply(_data);
  validateVector(_data, _compareDataResSum);
  precond.revert(_data);
  validateVector(_data, backup);
}

void PreconditionerTest::testValuePreconditioner ()
{
  preciceTrace("testValuePreconditioner()");
  std::vector<int> dims;
  dims.push_back(1);
  dims.push_back(2);
  dims.push_back(1);

  impl::ValuePreconditioner precond(dims);

  int numberOfRows = 8;
  precond.initialize(numberOfRows);
  DataValues backup = _data;

  preciceDebug("New iteration");
  //should change, since first timestep
  precond.update(false, _data, _res);
  validate(precond.requireNewQR());
  precond.newQRfulfilled();
  precond.apply(_data);
  validateVector(_data, _compareDataValue);
  precond.revert(_data);
  validateVector(_data, backup);

  //now no change
  preciceDebug("Another new iteration");
  precond.update(false, _data, _res);
  validate(not precond.requireNewQR());
  precond.apply(_data);
  validateVector(_data, _compareDataValue);
  precond.revert(_data);
  validateVector(_data, backup);

  preciceDebug("New timestep");
  //should change weights
  precond.update(true, _data*2, _res);
  validate(precond.requireNewQR());
  precond.newQRfulfilled();
}

void PreconditionerTest::testConstPreconditioner ()
{
  preciceTrace("testConstPreconditioner()");
  std::vector<int> dims;
  dims.push_back(1);
  dims.push_back(2);
  dims.push_back(1);

  std::vector<double> factors;
  factors.push_back(1e3);
  factors.push_back(2.0);
  factors.push_back(1e-5);

  impl::ConstantPreconditioner precond(dims, factors);

  int numberOfRows = 8;
  precond.initialize(numberOfRows); //new weights already computed here
  DataValues backup = _data;

  preciceDebug("New iteration");
  // should have no effect
  precond.update(false, _data, _res);
  validate(not precond.requireNewQR());
  precond.apply(_data);
  validateVector(_data, _compareDataConstant);
  precond.revert(_data);
  validateVector(_data, backup);

  preciceDebug("New timestep");
  //should not change weights
  precond.update(true, _data, _res);
  validate(not precond.requireNewQR());
  precond.apply(_data);
  validateVector(_data, _compareDataConstant);
  precond.revert(_data);
  validateVector(_data, backup);
}



void PreconditionerTest::testParallelMatrixScaling ()
{
#ifndef PRECICE_NO_MPI

  preciceTrace("testParallelMatrixScaling()");
  assertion ( utils::Parallel::getCommunicatorSize() == 4 );
  utils::Parallel::synchronizeProcesses();

  utils::MasterSlave::_communication = com::Communication::SharedPointer(new com::MPIDirectCommunication());

  utils::MasterSlave::_size = 4;
  utils::MasterSlave::_rank = utils::Parallel::getProcessRank();
  utils::MasterSlave::_slaveMode = false;
  utils::MasterSlave::_masterMode = false;


  //setup communication and master-slave layout
  if (utils::Parallel::getProcessRank() == 0){
    utils::Parallel::splitCommunicator( "PrecondMaster" );
    utils::MasterSlave::_masterMode = true;
    utils::MasterSlave::_communication->acceptConnection ( "PrecondMaster", "PrecondSlave", 0, 1);
    utils::MasterSlave::_communication->setRankOffset(1);
  }
  else if(utils::Parallel::getProcessRank() == 1){
    utils::Parallel::splitCommunicator( "PrecondSlave" );
    utils::MasterSlave::_slaveMode = true;
    utils::MasterSlave::_communication->requestConnection( "PrecondMaster", "PrecondSlave", 0, 3 );
  }
  else if(utils::Parallel::getProcessRank() == 2){
    utils::Parallel::splitCommunicator( "PrecondSlave");
    utils::MasterSlave::_slaveMode = true;
    utils::MasterSlave::_communication->requestConnection( "PrecondMaster", "PrecondSlave", 1, 3 );
  }
  else if(utils::Parallel::getProcessRank() == 3){
    utils::Parallel::splitCommunicator( "PrecondSlave");
    utils::MasterSlave::_slaveMode = true;
    utils::MasterSlave::_communication->requestConnection( "PrecondMaster", "PrecondSlave", 2, 3 );
  }



  //setup data
  int localN = -1;
  if (utils::Parallel::getProcessRank() == 0){
    localN = 2;
  }
  else if(utils::Parallel::getProcessRank() == 1){
    localN = 1;
  }
  else if(utils::Parallel::getProcessRank() == 2){
    localN = 0;
  }
  else if(utils::Parallel::getProcessRank() == 3){
    localN = 1;
  }

  int globalN = 4;

  Eigen::MatrixXd V(localN,2);
  Eigen::MatrixXd M(globalN,localN);
  Eigen::VectorXd x(localN);
  Eigen::MatrixXd V_back(localN,2);
  Eigen::MatrixXd M_back(globalN,localN);
  Eigen::VectorXd x_back(localN);


  if (utils::Parallel::getProcessRank() == 0){
    V(0,0) = 1.0;
    V(0,1) = 2.0;
    V(1,0) = 3.0;
    V(1,1) = 4.0;
    M(0,0) = 1.0;
    M(0,1) = 2.0;
    M(1,0) = 3.0;
    M(1,1) = 4.0;
    M(2,0) = 1.0;
    M(2,1) = 2.0;
    M(3,0) = 3.0;
    M(3,1) = 4.0;
    x(0) = 5.0;
    x(1) = 5.0;
  }
  else if(utils::Parallel::getProcessRank() == 1){
    V(0,0) = 5.0;
    V(0,1) = 6.0;
    M(0,0) = 1.0;
    M(1,0) = 2.0;
    M(2,0) = 3.0;
    M(3,0) = 4.0;
    x(0) = 5.0;
  }
  else if(utils::Parallel::getProcessRank() == 2){
  }
  else if(utils::Parallel::getProcessRank() == 3){
    V(0,0) = 7.0;
    V(0,1) = 8.0;
    M(0,0) = 1.0;
    M(1,0) = 2.0;
    M(2,0) = 3.0;
    M(3,0) = 4.0;
    x(0) = 5.0;
  }

  V_back = V;
  M_back = M;
  x_back = x;


  std::vector<int> dims;
  dims.push_back(1);


  impl::ValuePreconditioner precond(dims);
  precond.initialize(localN);
  precond.triggerGlobalWeights(globalN);
  precond.update(true,x,x);
  validate(precond.requireNewQR());

  precond.apply(V);

  for(int i=0; i<V.rows(); i++){
    for(int j=0; j<V.cols(); j++){
      validateWithParams2(tarch::la::equals(V(i,j), V_back(i,j)*0.1), V(i,j), V_back(i,j)*0.1);
    }
  }

  precond.revert(V);

  for(int i=0; i<V.rows(); i++){
    for(int j=0; j<V.cols(); j++){
      validateWithParams2(tarch::la::equals(V(i,j), V_back(i,j)), V(i,j), V_back(i,j));
    }
  }

  precond.apply(M, false);
  precond.revert(M, true);

  for(int i=0; i<M.rows(); i++){
    for(int j=0; j<M.cols(); j++){
      validateWithParams2(tarch::la::equals(M(i,j), M_back(i,j)), M(i,j), M_back(i,j));
    }
  }

  precond.revert(M, false);
  precond.apply(M, true);

  for(int i=0; i<M.rows(); i++){
    for(int j=0; j<M.cols(); j++){
      validateWithParams2(tarch::la::equals(M(i,j), M_back(i,j)), M(i,j), M_back(i,j));
    }
  }


  utils::MasterSlave::_slaveMode = false;
  utils::MasterSlave::_masterMode = false;
  utils::Parallel::clearGroups();
  utils::Parallel::synchronizeProcesses();
  utils::MasterSlave::_communication = nullptr;

#endif
}

void PreconditionerTest::validateVector (DataValues& data, DataValues& compare)
{
  validate(data.size()==compare.size());
  for(int i=0; i<data.size(); i++){
    validateWithParams2(tarch::la::equals(data(i), compare(i),1e-8), data(i), compare(i));
  }
}



}}} // namespace precice, cplscheme
