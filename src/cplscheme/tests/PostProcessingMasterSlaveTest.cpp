// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License

/*
 *  Created on: Aug 4, 2015
 *      Author: Klaudius Scheufele
 */

#ifndef PRECICE_NO_MPI

#include "PostProcessingMasterSlaveTest.hpp"
#include "utils/Globals.hpp"
#include "utils/Dimensions.hpp"
#include "utils/MasterSlave.hpp"
#include "mesh/Mesh.hpp"
#include "utils/Parallel.hpp"
#include "m2n/M2N.hpp"
#include "mapping/SharedPointer.hpp"
#include "cplscheme/impl/IQNILSPostProcessing.hpp"
#include "cplscheme/impl/MVQNPostProcessing.hpp"
#include "cplscheme/impl/BaseQNPostProcessing.hpp"
#include "cplscheme/ParallelCouplingScheme.hpp"
#include "cplscheme/impl/ConvergenceMeasure.hpp"
#include "cplscheme/impl/AbsoluteConvergenceMeasure.hpp"
#include "cplscheme/impl/MinIterationConvergenceMeasure.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "cplscheme/impl/SharedPointer.hpp"
#include "cplscheme/impl/ConstantPreconditioner.hpp"
#include "cplscheme/Constants.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "com/MPIPortsCommunication.hpp"
#include "tarch/la/Vector.h"
#include "tarch/la/WrappedVector.h"
#include <string.h>
#include "utils/EigenHelperFunctions.hpp"

#include "tarch/tests/TestCaseFactory.h"
registerTest(precice::cplscheme::tests::PostProcessingMasterSlaveTest)

namespace precice {
namespace cplscheme {
namespace tests {

using utils::Vector3D;

tarch::logging::Log PostProcessingMasterSlaveTest::
  _log ( "precice::cplscheme::tests::PostProcessingMasterSlaveTest" );

PostProcessingMasterSlaveTest:: PostProcessingMasterSlaveTest ()
:
  TestCase ( "cplscheme::tests::PostProcessingMasterSlaveTest" )
{}

void PostProcessingMasterSlaveTest:: run ()
{
  preciceTrace ( "run" );
  typedef utils::Parallel Par;
  if (Par::getCommunicatorSize() > 3){
    std::vector<int> ranksWanted;
    ranksWanted += 0, 1, 2 , 3;
    MPI_Comm comm = Par::getRestrictedCommunicator(ranksWanted);
    if (Par::getProcessRank() <= 3){
      Par::setGlobalCommunicator(comm);
      testMethod (testVIQNILSpp);
      Par::setGlobalCommunicator(Par::getCommunicatorWorld());
    }
    comm = Par::getRestrictedCommunicator(ranksWanted);
    if (Par::getProcessRank() <= 3){
      Par::setGlobalCommunicator(comm); //necessary to be able to re-initialize with different leading ranks
      // testMethod (testVIQNIMVJpp);
      Par::setGlobalCommunicator(Par::getCommunicatorWorld());
    }
  }
}

void PostProcessingMasterSlaveTest::testVIQNILSpp()
{

	preciceTrace ( "testVIQNPP" ); assertion ( utils::Parallel::getCommunicatorSize() == 4 );

	com::Communication::SharedPointer masterSlaveCom =
			com::Communication::SharedPointer(
					new com::MPIDirectCommunication());
	utils::MasterSlave::_communication = masterSlaveCom;

	utils::Parallel::synchronizeProcesses();

	if (utils::Parallel::getProcessRank() == 0) { //Master
		utils::Parallel::splitCommunicator("SOLIDZMaster");
	} else if (utils::Parallel::getProcessRank() == 1) { //Slave1
		utils::Parallel::splitCommunicator("SOLIDZSlaves");
	} else if (utils::Parallel::getProcessRank() == 2) { //Slave2
		utils::Parallel::splitCommunicator("SOLIDZSlaves");
	} else if (utils::Parallel::getProcessRank() == 3) { //Slave3
		utils::Parallel::splitCommunicator("SOLIDZSlaves");
	}

	if (utils::Parallel::getProcessRank() == 0) { //Master
		masterSlaveCom->acceptConnection("SOLIDZMaster", "SOLIDZSlaves", 0, 1);
		masterSlaveCom->setRankOffset(1);
	} else if (utils::Parallel::getProcessRank() == 1) { //Slave1
		masterSlaveCom->requestConnection("SOLIDZMaster", "SOLIDZSlaves", 0, 3);
	} else if (utils::Parallel::getProcessRank() == 2) { //Slave2
		masterSlaveCom->requestConnection("SOLIDZMaster", "SOLIDZSlaves", 1, 3);
	} else if (utils::Parallel::getProcessRank() == 3) { //Slave3
		masterSlaveCom->requestConnection("SOLIDZMaster", "SOLIDZSlaves", 2, 3);
	}


	double initialRelaxation = 0.01;
	int    maxIterationsUsed = 50;
	int    timestepsReused = 6;
	int filter = impl::BaseQNPostProcessing::QR1FILTER;
	double singularityLimit = 1e-10;
	bool enforceInitialRelaxation = false;
	std::vector<int> dataIDs;
	dataIDs.push_back(0);
	dataIDs.push_back(1);
	std::vector<double> factors;
  factors.resize(2,1.0);
  std::vector<int> dims;
  dims.resize(2,1);
  impl::PtrPreconditioner prec(new impl::ConstantPreconditioner(dims,factors));
	std::vector<int> vertexOffsets {4, 8, 8 , 10};

	mesh::PtrMesh dummyMesh ( new mesh::Mesh("dummyMesh", 3, false) );
	dummyMesh->setVertexOffsets(vertexOffsets);

	cplscheme::impl::IQNILSPostProcessing pp(initialRelaxation, enforceInitialRelaxation, maxIterationsUsed,
										   timestepsReused, filter, singularityLimit, dataIDs, prec);

	Eigen::VectorXd dvalues;
	Eigen::VectorXd dcol1;
	Eigen::VectorXd fvalues;
	Eigen::VectorXd fcol1;

	DataMap data;


	if (utils::Parallel::getProcessRank() == 0) { //Master
		utils::MasterSlave::_rank = 0;
		utils::MasterSlave::_size = 4;
		utils::MasterSlave::_slaveMode = false;
		utils::MasterSlave::_masterMode = true;

		/**
		 * processor with 4 vertices
		 */

		//init displacements
		Eigen::VectorXd insert(4); insert << 1.0, 2.0, 3.0, 4.0;
    utils::append(dvalues, insert);
    insert << 1.0, 1.0, 1.0, 1.0;
    utils::append(dcol1, insert);
		//dvalues.append(1.0); dvalues.append(2.0); dvalues.append(3.0); dvalues.append(4.0);
		//dcol1.append(1.0); dcol1.append(1.0); dcol1.append(1.0); dcol1.append(1.0);

		PtrCouplingData dpcd(new CouplingData(&dvalues,dummyMesh,false,1));

		//init forces
    insert << 0.1, 0.1, 0.1, 0.1;
    utils::append(fvalues, insert);
    insert << 0.2, 0.2, 0.2, 0.2;
    utils::append(fcol1, insert);
		//fvalues.append(0.1); fvalues.append(0.1); fvalues.append(0.1); fvalues.append(0.1);
		//fcol1.append(0.2); fcol1.append(0.2); fcol1.append(0.2); fcol1.append(0.2);

		PtrCouplingData fpcd(new CouplingData(&fvalues,dummyMesh,false,1));

		data.insert(std::pair<int,PtrCouplingData>(0,dpcd));
		data.insert(std::pair<int,PtrCouplingData>(1,fpcd));

		pp.initialize(data);

		dpcd->oldValues.col(0) = dcol1;
		fpcd->oldValues.col(0) = fcol1;

	} else if (utils::Parallel::getProcessRank() == 1) { //Slave1
		utils::MasterSlave::_rank = 1;
		utils::MasterSlave::_size = 4;
		utils::MasterSlave::_slaveMode = true;
		utils::MasterSlave::_masterMode = false;

		/**
		 * processor with 4 vertices
		 */

		//init displacements
		Eigen::VectorXd insert(4); insert << 5.0, 6.0, 7.0, 8.0;
		utils::append(dvalues, insert);
		insert << 1.0, 1.0, 1.0, 1.0;
		utils::append(dcol1, insert);

		//dvalues.append(5.0); dvalues.append(6.0); dvalues.append(7.0); dvalues.append(8.0);
		//dcol1.append(1.0); dcol1.append(1.0); dcol1.append(1.0); dcol1.append(1.0);

		PtrCouplingData dpcd(new CouplingData(&dvalues,dummyMesh,false,1));

		//init forces
		insert << 0.1, 0.1, 0.1, 0.1;
    utils::append(fvalues, insert);
    insert << 0.2, 0.2, 0.2, 0.2;
    utils::append(fcol1, insert);

    //fvalues.append(0.1); fvalues.append(0.1); fvalues.append(0.1); fvalues.append(0.1);
		//fcol1.append(0.2); fcol1.append(0.2); fcol1.append(0.2); fcol1.append(0.2);

		PtrCouplingData fpcd(new CouplingData(&fvalues,dummyMesh,false,1));

		data.insert(std::pair<int,PtrCouplingData>(0,dpcd));
		data.insert(std::pair<int,PtrCouplingData>(1,fpcd));

		pp.initialize(data);

		dpcd->oldValues.col(0) = dcol1;
		fpcd->oldValues.col(0) = fcol1;

	} else if (utils::Parallel::getProcessRank() == 2) { //Slave2
		utils::MasterSlave::_rank = 2;
		utils::MasterSlave::_size = 4;
		utils::MasterSlave::_slaveMode = true;
		utils::MasterSlave::_masterMode = false;

		/**
		 * processor with no vertices
		 */

		//init displacements
		PtrCouplingData dpcd(new CouplingData(&dvalues,dummyMesh,false,1));

		//init forces
		PtrCouplingData fpcd(new CouplingData(&fvalues,dummyMesh,false,1));

		data.insert(std::pair<int,PtrCouplingData>(0,dpcd));
		data.insert(std::pair<int,PtrCouplingData>(1,fpcd));

		pp.initialize(data);

		dpcd->oldValues.col(0) = dcol1;
		fpcd->oldValues.col(0) = fcol1;

	} else if (utils::Parallel::getProcessRank() == 3) { //Slave3
		utils::MasterSlave::_rank = 3;
		utils::MasterSlave::_size = 4;
		utils::MasterSlave::_slaveMode = true;
		utils::MasterSlave::_masterMode = false;

		/**
		 * processor with 2 vertices
		 */

		//init displacements
    Eigen::VectorXd insert(2); insert << 1.0, 2.0;
    utils::append(dvalues, insert);
    insert << 1.0, 1.0;
    utils::append(dcol1, insert);
		//dvalues.append(1.0); dvalues.append(2.0);
		//dcol1.append(1.0); dcol1.append(1.0);

		PtrCouplingData dpcd(new CouplingData(&dvalues,dummyMesh,false,1));

		//init forces
    insert << 0.1, 0.1;
    utils::append(fvalues, insert);
    insert << 0.2, 0.2;
    utils::append(fcol1, insert);
		//fvalues.append(0.1); fvalues.append(0.1);
		//fcol1.append(0.2); fcol1.append(0.2);

		PtrCouplingData fpcd(new CouplingData(&fvalues,dummyMesh,false,1));

		data.insert(std::pair<int,PtrCouplingData>(0,dpcd));
		data.insert(std::pair<int,PtrCouplingData>(1,fpcd));

		pp.initialize(data);

		dpcd->oldValues.col(0) = dcol1;
		fpcd->oldValues.col(0) = fcol1;
	}

	utils::Parallel::synchronizeProcesses();
	pp.performPostProcessing(data);
	utils::Parallel::synchronizeProcesses();

	Eigen::VectorXd newdvalues;
	if (utils::Parallel::getProcessRank() == 0) { //Master

		validateWithParams1(tarch::la::equals((*data.at(0)->values)(0), 1.00), (*data.at(0)->values)(0));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(1), 1.01), (*data.at(0)->values)(1));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(2), 1.02), (*data.at(0)->values)(2));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(3), 1.03), (*data.at(0)->values)(3));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(0), 0.199), (*data.at(1)->values)(0));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(1), 0.199), (*data.at(1)->values)(1));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(2), 0.199), (*data.at(1)->values)(2));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(3), 0.199), (*data.at(1)->values)(3));

		/*
		std::cout<<"  Master:"<<std::endl;
		std::cout<<"Master: (*data.at(0)->values)(0): "<<(*data.at(0)->values)(0)<<std::endl;
		std::cout<<"Master: (*data.at(0)->values)(1): "<<(*data.at(0)->values)(1)<<std::endl;
		std::cout<<"Master: (*data.at(0)->values)(2): "<<(*data.at(0)->values)(2)<<std::endl;
		std::cout<<"Master: (*data.at(0)->values)(3): "<<(*data.at(0)->values)(3)<<std::endl;
		std::cout<<"Master: (*data.at(1)->values)(0): "<<(*data.at(1)->values)(0)<<std::endl;
		std::cout<<"Master: (*data.at(1)->values)(1): "<<(*data.at(1)->values)(1)<<std::endl;
		std::cout<<"Master: (*data.at(1)->values)(2): "<<(*data.at(1)->values)(2)<<std::endl;
		std::cout<<"Master: (*data.at(1)->values)(3): "<<(*data.at(1)->values)(3)<<std::endl;
		*/

		utils::append(newdvalues, 10.0); utils::append(newdvalues, 10.0); utils::append(newdvalues, 10.0); utils::append(newdvalues, 10.0);
//		newdvalues.append(10.0); newdvalues.append(10.0); newdvalues.append(10.0); newdvalues.append(10.0);

	} else if (utils::Parallel::getProcessRank() == 1) { //Slave1

		validateWithParams1(tarch::la::equals((*data.at(0)->values)(0), 1.04), (*data.at(0)->values)(0));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(1), 1.05), (*data.at(0)->values)(1));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(2), 1.06), (*data.at(0)->values)(2));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(3), 1.07), (*data.at(0)->values)(3));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(0), 0.199), (*data.at(1)->values)(0));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(1), 0.199), (*data.at(1)->values)(1));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(2), 0.199), (*data.at(1)->values)(2));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(3), 0.199), (*data.at(1)->values)(3));

		/*
		std::cout<<"  Slave 1:"<<std::endl;
		std::cout<<"Slave 1: (*data.at(0)->values)(0): "<<(*data.at(0)->values)(0)<<std::endl;
		std::cout<<"Slave 1: (*data.at(0)->values)(1): "<<(*data.at(0)->values)(1)<<std::endl;
		std::cout<<"Slave 1: (*data.at(0)->values)(2): "<<(*data.at(0)->values)(2)<<std::endl;
		std::cout<<"Slave 1: (*data.at(0)->values)(3): "<<(*data.at(0)->values)(3)<<std::endl;
		std::cout<<"Slave 1: (*data.at(1)->values)(0): "<<(*data.at(1)->values)(0)<<std::endl;
		std::cout<<"Slave 1: (*data.at(1)->values)(1): "<<(*data.at(1)->values)(1)<<std::endl;
		std::cout<<"Slave 1: (*data.at(1)->values)(2): "<<(*data.at(1)->values)(2)<<std::endl;
		std::cout<<"Slave 1: (*data.at(1)->values)(3): "<<(*data.at(1)->values)(3)<<std::endl;
		 */

		utils::append(newdvalues, 10.0); utils::append(newdvalues, 10.0); utils::append(newdvalues, 10.0); utils::append(newdvalues, 10.0);
//		newdvalues.append(10.0); newdvalues.append(10.0); newdvalues.append(10.0); newdvalues.append(10.0);

	} else if (utils::Parallel::getProcessRank() == 2) { //Slave2
		// empty proc

	} else if (utils::Parallel::getProcessRank() == 3) { //Slave3

		validateWithParams1(tarch::la::equals((*data.at(0)->values)(0), 1.00), (*data.at(0)->values)(0));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(1), 1.01), (*data.at(0)->values)(1));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(0), 0.199), (*data.at(1)->values)(0));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(1), 0.199), (*data.at(1)->values)(1));

		/*
		std::cout<<"  Slave 3:"<<std::endl;
		std::cout<<"Slave 3: (*data.at(0)->values)(0): "<<(*data.at(0)->values)(0)<<std::endl;
		std::cout<<"Slave 3: (*data.at(0)->values)(1): "<<(*data.at(0)->values)(1)<<std::endl;
		std::cout<<"Slave 3: (*data.at(1)->values)(0): "<<(*data.at(1)->values)(0)<<std::endl;
		std::cout<<"Slave 3: (*data.at(1)->values)(1): "<<(*data.at(1)->values)(1)<<std::endl;
		 */

		utils::append(newdvalues, 10.0); utils::append(newdvalues, 10.0);
		//newdvalues.append(10.0); newdvalues.append(10.0);
	}

	data.begin()->second->values = &newdvalues;

	utils::Parallel::synchronizeProcesses();
	pp.performPostProcessing(data);
	utils::Parallel::synchronizeProcesses();


	if (utils::Parallel::getProcessRank() == 0) { //Master

		validateWithParams1(tarch::la::equals((*data.at(0)->values)(0), -1.51483105223442748866e+00), (*data.at(0)->values)(0));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(1), -2.35405379763935940218e-01), (*data.at(0)->values)(1));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(2), 1.04402029270655560822e+00), (*data.at(0)->values)(2));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(3), 2.32344596517704804484e+00), (*data.at(0)->values)(3));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(0), 7.23368584254212854123e-02), (*data.at(1)->values)(0));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(1), 7.23368584254212854123e-02), (*data.at(1)->values)(1));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(2), 7.23368584254212854123e-02), (*data.at(1)->values)(2));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(3), 7.23368584254212854123e-02), (*data.at(1)->values)(3));

/*
		std::cout<<"  Master:"<<std::endl;
		std::cout<<"Master: (*data.at(0)->values)(0): "<<(*data.at(0)->values)(0)<<std::endl;
		std::cout<<"Master: (*data.at(0)->values)(1): "<<(*data.at(0)->values)(1)<<std::endl;
		std::cout<<"Master: (*data.at(0)->values)(2): "<<(*data.at(0)->values)(2)<<std::endl;
		std::cout<<"Master: (*data.at(0)->values)(3): "<<(*data.at(0)->values)(3)<<std::endl;
		std::cout<<"Master: (*data.at(1)->values)(0): "<<(*data.at(1)->values)(0)<<std::endl;
		std::cout<<"Master: (*data.at(1)->values)(1): "<<(*data.at(1)->values)(1)<<std::endl;
		std::cout<<"Master: (*data.at(1)->values)(2): "<<(*data.at(1)->values)(2)<<std::endl;
		std::cout<<"Master: (*data.at(1)->values)(3): "<<(*data.at(1)->values)(3)<<std::endl;
*/

	} else if (utils::Parallel::getProcessRank() == 1) { //Slave1
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(0), 3.60287163764754048145e+00), (*data.at(0)->values)(0));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(1), 4.88229731011803202989e+00), (*data.at(0)->values)(1));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(2), 6.16172298258852357833e+00), (*data.at(0)->values)(2));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(3), 7.44114865505901601495e+00), (*data.at(0)->values)(3));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(0), 7.23368584254212854123e-02), (*data.at(1)->values)(0));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(1), 7.23368584254212854123e-02), (*data.at(1)->values)(1));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(2), 7.23368584254212854123e-02), (*data.at(1)->values)(2));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(3), 7.23368584254212854123e-02), (*data.at(1)->values)(3));

		/*
		std::cout<<"  Slave 1:"<<std::endl;
		std::cout<<"Slave 1: (*data.at(0)->values)(0): "<<(*data.at(0)->values)(0)<<std::endl;
		std::cout<<"Slave 1: (*data.at(0)->values)(1): "<<(*data.at(0)->values)(1)<<std::endl;
		std::cout<<"Slave 1: (*data.at(0)->values)(2): "<<(*data.at(0)->values)(2)<<std::endl;
		std::cout<<"Slave 1: (*data.at(0)->values)(3): "<<(*data.at(0)->values)(3)<<std::endl;
		std::cout<<"Slave 1: (*data.at(1)->values)(0): "<<(*data.at(1)->values)(0)<<std::endl;
		std::cout<<"Slave 1: (*data.at(1)->values)(1): "<<(*data.at(1)->values)(1)<<std::endl;
		std::cout<<"Slave 1: (*data.at(1)->values)(2): "<<(*data.at(1)->values)(2)<<std::endl;
		std::cout<<"Slave 1: (*data.at(1)->values)(3): "<<(*data.at(1)->values)(3)<<std::endl;
		*/
	} else if (utils::Parallel::getProcessRank() == 2) { //Slave2
		// empty proc

	} else if (utils::Parallel::getProcessRank() == 3) { //Slave3

		validateWithParams1(tarch::la::equals((*data.at(0)->values)(0), -1.51483105223442748866e+00), (*data.at(0)->values)(0));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(1), -2.35405379763935940218e-01), (*data.at(0)->values)(1));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(0), 7.23368584254212854123e-02), (*data.at(1)->values)(0));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(1), 7.23368584254212854123e-02), (*data.at(1)->values)(1));

		/*
		std::cout<<"  Slave 3:"<<std::endl;
		std::cout<<"Slave 3: (*data.at(0)->values)(0): "<<(*data.at(0)->values)(0)<<std::endl;
		std::cout<<"Slave 3: (*data.at(0)->values)(1): "<<(*data.at(0)->values)(1)<<std::endl;
		std::cout<<"Slave 3: (*data.at(1)->values)(0): "<<(*data.at(1)->values)(0)<<std::endl;
		std::cout<<"Slave 3: (*data.at(1)->values)(1): "<<(*data.at(1)->values)(1)<<std::endl;
		*/
	}


	utils::Parallel::synchronizeProcesses();
	utils::MasterSlave::_slaveMode = false;
	utils::MasterSlave::_masterMode = false;
	utils::Parallel::clearGroups();
}


void PostProcessingMasterSlaveTest::testVIQNIMVJpp()
{
	preciceTrace ( "testVIQNIMVJpp" ); assertion ( utils::Parallel::getCommunicatorSize() == 4 );

	com::Communication::SharedPointer masterSlaveCom = com::Communication::SharedPointer(new com::MPIPortsCommunication("."));
	utils::MasterSlave::_communication = masterSlaveCom;

	utils::Parallel::synchronizeProcesses();

	if (utils::Parallel::getProcessRank() == 0) { //Master
		utils::Parallel::splitCommunicator("SOLIDZMaster");
	} else if (utils::Parallel::getProcessRank() == 1) { //Slave1
		utils::Parallel::splitCommunicator("SOLIDZSlaves");
	} else if (utils::Parallel::getProcessRank() == 2) { //Slave2
		utils::Parallel::splitCommunicator("SOLIDZSlaves");
	} else if (utils::Parallel::getProcessRank() == 3) { //Slave3
		utils::Parallel::splitCommunicator("SOLIDZSlaves");
	}


	if (utils::Parallel::getProcessRank() == 0) { //Master
		masterSlaveCom->acceptConnection("SOLIDZMaster", "SOLIDZSlaves", 0, 1);
		masterSlaveCom->setRankOffset(1);
	} else if (utils::Parallel::getProcessRank() == 1) { //Slave1
		masterSlaveCom->requestConnection("SOLIDZMaster", "SOLIDZSlaves", 0, 3);
	} else if (utils::Parallel::getProcessRank() == 2) { //Slave2
		masterSlaveCom->requestConnection("SOLIDZMaster", "SOLIDZSlaves", 1, 3);
	} else if (utils::Parallel::getProcessRank() == 3) { //Slave3
		masterSlaveCom->requestConnection("SOLIDZMaster", "SOLIDZSlaves", 2, 3);
	}


	double initialRelaxation = 0.01;
	int    maxIterationsUsed = 50;
	int    timestepsReused = 6;
	int filter = impl::BaseQNPostProcessing::QR1FILTER;
	int restartType = impl::MVQNPostProcessing::NO_RESTART;
	int chunkSize = 0;
	int reusedTimeStepsAtRestart = 0;
	double singularityLimit = 1e-10;
	double svdTruncationEps = 0.0;
	bool enforceInitialRelaxation = false;
	bool alwaysBuildJacobian = false;

	std::vector<int> dataIDs;
	dataIDs.push_back(0);
	dataIDs.push_back(1);
	std::vector<double> factors;
  factors.resize(2,1.0);
  std::vector<int> dims;
  dims.resize(2,1);
  impl::PtrPreconditioner prec(new impl::ConstantPreconditioner(dims,factors));
	std::vector<int> vertexOffsets {4, 8, 8 , 10};

	mesh::PtrMesh dummyMesh ( new mesh::Mesh("dummyMesh", 3, false) );
	dummyMesh->setVertexOffsets(vertexOffsets);

	cplscheme::impl::MVQNPostProcessing pp(initialRelaxation, enforceInitialRelaxation, maxIterationsUsed,
									   timestepsReused, filter, singularityLimit, dataIDs, prec, alwaysBuildJacobian,
									   restartType, chunkSize, reusedTimeStepsAtRestart, svdTruncationEps);

	Eigen::VectorXd dvalues;
	Eigen::VectorXd dcol1;
	Eigen::VectorXd fvalues;
	Eigen::VectorXd fcol1;

	DataMap data;


	if (utils::Parallel::getProcessRank() == 0) { //Master
		utils::MasterSlave::_rank = 0;
		utils::MasterSlave::_size = 4;
		utils::MasterSlave::_slaveMode = false;
		utils::MasterSlave::_masterMode = true;

		/**
		 * processor with 4 vertices
		 */

		//init displacements
    Eigen::VectorXd insert(4); insert << 1.0, 2.0, 3.0, 4.0;
    utils::append(dvalues, insert);
    insert << 1.0, 1.0, 1.0, 1.0;
    utils::append(dcol1, insert);
    //dvalues.append(1.0); dvalues.append(2.0); dvalues.append(3.0); dvalues.append(4.0);
    //dcol1.append(1.0); dcol1.append(1.0); dcol1.append(1.0); dcol1.append(1.0);

    PtrCouplingData dpcd(new CouplingData(&dvalues,dummyMesh,false,1));

    //init forces
    insert << 0.1, 0.1, 0.1, 0.1;
    utils::append(fvalues, insert);
    insert << 0.2, 0.2, 0.2, 0.2;
    utils::append(fcol1, insert);
    //fvalues.append(0.1); fvalues.append(0.1); fvalues.append(0.1); fvalues.append(0.1);
    //fcol1.append(0.2); fcol1.append(0.2); fcol1.append(0.2); fcol1.append(0.2);

		PtrCouplingData fpcd(new CouplingData(&fvalues,dummyMesh,false,1));

		data.insert(std::pair<int,PtrCouplingData>(0,dpcd));
		data.insert(std::pair<int,PtrCouplingData>(1,fpcd));

		pp.initialize(data);

		dpcd->oldValues.col(0) = dcol1;
		fpcd->oldValues.col(0) = fcol1;

	} else if (utils::Parallel::getProcessRank() == 1) { //Slave1
		utils::MasterSlave::_rank = 1;
		utils::MasterSlave::_size = 4;
		utils::MasterSlave::_slaveMode = true;
		utils::MasterSlave::_masterMode = false;

		/**
		 * processor with 4 vertices
		 */

		//init displacements
		Eigen::VectorXd insert(4); insert << 5.0, 6.0, 7.0, 8.0;
    utils::append(dvalues, insert);
    insert << 1.0, 1.0, 1.0, 1.0;
    utils::append(dcol1, insert);
		//dvalues.append(5.0); dvalues.append(6.0); dvalues.append(7.0); dvalues.append(8.0);
		//dcol1.append(1.0); dcol1.append(1.0); dcol1.append(1.0); dcol1.append(1.0);

		PtrCouplingData dpcd(new CouplingData(&dvalues,dummyMesh,false,1));

		//init forces
		insert << 0.1, 0.1, 0.1, 0.1;
    utils::append(fvalues, insert);
    insert << 0.2, 0.2, 0.2, 0.2;
    utils::append(fcol1, insert);
		//fvalues.append(0.1); fvalues.append(0.1); fvalues.append(0.1); fvalues.append(0.1);
		//fcol1.append(0.2); fcol1.append(0.2); fcol1.append(0.2); fcol1.append(0.2);

		PtrCouplingData fpcd(new CouplingData(&fvalues,dummyMesh,false,1));

		data.insert(std::pair<int,PtrCouplingData>(0,dpcd));
		data.insert(std::pair<int,PtrCouplingData>(1,fpcd));

		pp.initialize(data);

		dpcd->oldValues.col(0) = dcol1;
		fpcd->oldValues.col(0) = fcol1;

	} else if (utils::Parallel::getProcessRank() == 2) { //Slave2
		utils::MasterSlave::_rank = 2;
		utils::MasterSlave::_size = 4;
		utils::MasterSlave::_slaveMode = true;
		utils::MasterSlave::_masterMode = false;

		/**
		 * processor with no vertices
		 */

		//init displacements
		PtrCouplingData dpcd(new CouplingData(&dvalues,dummyMesh,false,1));

		//init forces
		PtrCouplingData fpcd(new CouplingData(&fvalues,dummyMesh,false,1));

		data.insert(std::pair<int,PtrCouplingData>(0,dpcd));
		data.insert(std::pair<int,PtrCouplingData>(1,fpcd));

		pp.initialize(data);

		dpcd->oldValues.col(0) = dcol1;
		fpcd->oldValues.col(0) = fcol1;

	} else if (utils::Parallel::getProcessRank() == 3) { //Slave3
		utils::MasterSlave::_rank = 3;
		utils::MasterSlave::_size = 4;
		utils::MasterSlave::_slaveMode = true;
		utils::MasterSlave::_masterMode = false;

		/**
		 * processor with 2 vertices
		 */

		//init displacements
		Eigen::VectorXd insert(2); insert << 1.0, 2.0;
    utils::append(dvalues, insert);
    insert << 1.0, 1.0;
    utils::append(dcol1, insert);
		//dvalues.append(1.0); dvalues.append(2.0);
		//dcol1.append(1.0); dcol1.append(1.0);

		PtrCouplingData dpcd(new CouplingData(&dvalues,dummyMesh,false,1));

		//init forces
		insert << 0.1, 0.1;
    utils::append(fvalues, insert);
    insert << 0.2, 0.2;
    utils::append(fcol1, insert);
		//fvalues.append(0.1); fvalues.append(0.1);
		//fcol1.append(0.2); fcol1.append(0.2);

		PtrCouplingData fpcd(new CouplingData(&fvalues,dummyMesh,false,1));

		data.insert(std::pair<int,PtrCouplingData>(0,dpcd));
		data.insert(std::pair<int,PtrCouplingData>(1,fpcd));

		pp.initialize(data);

		dpcd->oldValues.col(0) = dcol1;
		fpcd->oldValues.col(0) = fcol1;
	}

	pp.performPostProcessing(data);

	Eigen::VectorXd newdvalues;
	if (utils::Parallel::getProcessRank() == 0) { //Master

		validateWithParams1(tarch::la::equals((*data.at(0)->values)(0), 1.00000000000000000000e+00), (*data.at(0)->values)(0));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(1), 1.01000000000000000888e+00), (*data.at(0)->values)(1));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(2), 1.02000000000000001776e+00), (*data.at(0)->values)(2));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(3), 1.03000000000000002665e+00), (*data.at(0)->values)(3));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(0), 1.99000000000000010214e-01), (*data.at(1)->values)(0));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(1), 1.99000000000000010214e-01), (*data.at(1)->values)(1));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(2), 1.99000000000000010214e-01), (*data.at(1)->values)(2));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(3), 1.99000000000000010214e-01), (*data.at(1)->values)(3));
		/*
		std::cout<<"  Master:"<<std::endl;
		std::cout<<"Master: (*data.at(0)->values)(0): "<<(*data.at(0)->values)(0)<<std::endl;
		std::cout<<"Master: (*data.at(0)->values)(1): "<<(*data.at(0)->values)(1)<<std::endl;
		std::cout<<"Master: (*data.at(0)->values)(2): "<<(*data.at(0)->values)(2)<<std::endl;
		std::cout<<"Master: (*data.at(0)->values)(3): "<<(*data.at(0)->values)(3)<<std::endl;
		std::cout<<"Master: (*data.at(1)->values)(0): "<<(*data.at(1)->values)(0)<<std::endl;
		std::cout<<"Master: (*data.at(1)->values)(1): "<<(*data.at(1)->values)(1)<<std::endl;
		std::cout<<"Master: (*data.at(1)->values)(2): "<<(*data.at(1)->values)(2)<<std::endl;
		std::cout<<"Master: (*data.at(1)->values)(3): "<<(*data.at(1)->values)(3)<<std::endl;
		*/

		utils::append(newdvalues, 10.0); utils::append(newdvalues, 10.0); utils::append(newdvalues, 10.0); utils::append(newdvalues, 10.0);
//		newdvalues.append(10.0); newdvalues.append(10.0); newdvalues.append(10.0); newdvalues.append(10.0);

	} else if (utils::Parallel::getProcessRank() == 1) { //Slave1

		validateWithParams1(tarch::la::equals((*data.at(0)->values)(0), 1.04000000000000003553e+00), (*data.at(0)->values)(0));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(1), 1.05000000000000004441e+00), (*data.at(0)->values)(1));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(2), 1.06000000000000005329e+00), (*data.at(0)->values)(2));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(3), 1.07000000000000006217e+00), (*data.at(0)->values)(3));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(0), 1.99000000000000010214e-01), (*data.at(1)->values)(0));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(1), 1.99000000000000010214e-01), (*data.at(1)->values)(1));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(2), 1.99000000000000010214e-01), (*data.at(1)->values)(2));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(3), 1.99000000000000010214e-01), (*data.at(1)->values)(3));
		/*
		std::cout<<"  Slave 1:"<<std::endl;
		std::cout<<"Slave 1: (*data.at(0)->values)(0): "<<(*data.at(0)->values)(0)<<std::endl;
		std::cout<<"Slave 1: (*data.at(0)->values)(1): "<<(*data.at(0)->values)(1)<<std::endl;
		std::cout<<"Slave 1: (*data.at(0)->values)(2): "<<(*data.at(0)->values)(2)<<std::endl;
		std::cout<<"Slave 1: (*data.at(0)->values)(3): "<<(*data.at(0)->values)(3)<<std::endl;
		std::cout<<"Slave 1: (*data.at(1)->values)(0): "<<(*data.at(1)->values)(0)<<std::endl;
		std::cout<<"Slave 1: (*data.at(1)->values)(1): "<<(*data.at(1)->values)(1)<<std::endl;
		std::cout<<"Slave 1: (*data.at(1)->values)(2): "<<(*data.at(1)->values)(2)<<std::endl;
		std::cout<<"Slave 1: (*data.at(1)->values)(3): "<<(*data.at(1)->values)(3)<<std::endl;
		*/

		utils::append(newdvalues, 10.0); utils::append(newdvalues, 10.0); utils::append(newdvalues, 10.0); utils::append(newdvalues, 10.0);
//		newdvalues.append(10.0); newdvalues.append(10.0); newdvalues.append(10.0); newdvalues.append(10.0);

	} else if (utils::Parallel::getProcessRank() == 2) { //Slave2
		// empty proc

	} else if (utils::Parallel::getProcessRank() == 3) { //Slave3

		validateWithParams1(tarch::la::equals((*data.at(0)->values)(0), 1.00000000000000000000e+00), (*data.at(0)->values)(0));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(1), 1.01000000000000000888e+00), (*data.at(0)->values)(1));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(0), 1.99000000000000010214e-01), (*data.at(1)->values)(0));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(1), 1.99000000000000010214e-01), (*data.at(1)->values)(1));
		/*
		std::cout<<"  Slave 3:"<<std::endl;
		std::cout<<"Slave 3: (*data.at(0)->values)(0): "<<(*data.at(0)->values)(0)<<std::endl;
		std::cout<<"Slave 3: (*data.at(0)->values)(1): "<<(*data.at(0)->values)(1)<<std::endl;
		std::cout<<"Slave 3: (*data.at(1)->values)(0): "<<(*data.at(1)->values)(0)<<std::endl;
		std::cout<<"Slave 3: (*data.at(1)->values)(1): "<<(*data.at(1)->values)(1)<<std::endl;
		*/

		utils::append(newdvalues, 10.0); utils::append(newdvalues, 10.0);
//		newdvalues.append(10.0); newdvalues.append(10.0);
	}

	data.begin()->second->values = &newdvalues;
	pp.performPostProcessing(data);

	if (utils::Parallel::getProcessRank() == 0) { //Master
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(0), -1.51483105223442748866e+00), (*data.at(0)->values)(0));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(1), -2.35405379763935940218e-01), (*data.at(0)->values)(1));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(2), 1.04402029270655738458e+00), (*data.at(0)->values)(2));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(3), 2.32344596517704893301e+00), (*data.at(0)->values)(3));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(0), 7.23368584254213131679e-02), (*data.at(1)->values)(0));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(1), 7.23368584254213131679e-02), (*data.at(1)->values)(1));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(2), 7.23368584254213131679e-02), (*data.at(1)->values)(2));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(3), 7.23368584254213131679e-02), (*data.at(1)->values)(3));
/*
		std::cout<<"  Master:"<<std::endl;
		std::cout<<"Master: (*data.at(0)->values)(0): "<<(*data.at(0)->values)(0)<<std::endl;
		std::cout<<"Master: (*data.at(0)->values)(1): "<<(*data.at(0)->values)(1)<<std::endl;
		std::cout<<"Master: (*data.at(0)->values)(2): "<<(*data.at(0)->values)(2)<<std::endl;
		std::cout<<"Master: (*data.at(0)->values)(3): "<<(*data.at(0)->values)(3)<<std::endl;
		std::cout<<"Master: (*data.at(1)->values)(0): "<<(*data.at(1)->values)(0)<<std::endl;
		std::cout<<"Master: (*data.at(1)->values)(1): "<<(*data.at(1)->values)(1)<<std::endl;
		std::cout<<"Master: (*data.at(1)->values)(2): "<<(*data.at(1)->values)(2)<<std::endl;
		std::cout<<"Master: (*data.at(1)->values)(3): "<<(*data.at(1)->values)(3)<<std::endl;
*/

	} else if (utils::Parallel::getProcessRank() == 1) { //Slave1
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(0), 3.60287163764754048145e+00), (*data.at(0)->values)(0));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(1), 4.88229731011803202989e+00), (*data.at(0)->values)(1));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(2), 6.16172298258852446651e+00), (*data.at(0)->values)(2));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(3), 7.44114865505901601495e+00), (*data.at(0)->values)(3));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(0), 7.23368584254213131679e-02), (*data.at(1)->values)(0));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(1), 7.23368584254213131679e-02), (*data.at(1)->values)(1));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(2), 7.23368584254213131679e-02), (*data.at(1)->values)(2));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(3), 7.23368584254213131679e-02), (*data.at(1)->values)(3));
		/*
		std::cout<<"  Slave 1:"<<std::endl;
		std::cout<<"Slave 1: (*data.at(0)->values)(0): "<<(*data.at(0)->values)(0)<<std::endl;
		std::cout<<"Slave 1: (*data.at(0)->values)(1): "<<(*data.at(0)->values)(1)<<std::endl;
		std::cout<<"Slave 1: (*data.at(0)->values)(2): "<<(*data.at(0)->values)(2)<<std::endl;
		std::cout<<"Slave 1: (*data.at(0)->values)(3): "<<(*data.at(0)->values)(3)<<std::endl;
		std::cout<<"Slave 1: (*data.at(1)->values)(0): "<<(*data.at(1)->values)(0)<<std::endl;
		std::cout<<"Slave 1: (*data.at(1)->values)(1): "<<(*data.at(1)->values)(1)<<std::endl;
		std::cout<<"Slave 1: (*data.at(1)->values)(2): "<<(*data.at(1)->values)(2)<<std::endl;
		std::cout<<"Slave 1: (*data.at(1)->values)(3): "<<(*data.at(1)->values)(3)<<std::endl;
		 */
	} else if (utils::Parallel::getProcessRank() == 2) { //Slave2
		// empty proc

	} else if (utils::Parallel::getProcessRank() == 3) { //Slave3

		validateWithParams1(tarch::la::equals((*data.at(0)->values)(0), -1.51483105223442748866e+00), (*data.at(0)->values)(0));
		validateWithParams1(tarch::la::equals((*data.at(0)->values)(1), -2.35405379763935940218e-01), (*data.at(0)->values)(1));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(0), 7.23368584254213131679e-02), (*data.at(1)->values)(0));
		validateWithParams1(tarch::la::equals((*data.at(1)->values)(1), 7.23368584254213131679e-02), (*data.at(1)->values)(1));
		/*
		std::cout<<"  Slave 3:"<<std::endl;
		std::cout<<"Slave 3: (*data.at(0)->values)(0): "<<(*data.at(0)->values)(0)<<std::endl;
		std::cout<<"Slave 3: (*data.at(0)->values)(1): "<<(*data.at(0)->values)(1)<<std::endl;
		std::cout<<"Slave 3: (*data.at(1)->values)(0): "<<(*data.at(1)->values)(0)<<std::endl;
		std::cout<<"Slave 3: (*data.at(1)->values)(1): "<<(*data.at(1)->values)(1)<<std::endl;
		*/
	}

	utils::MasterSlave::_communication->closeConnection();
	utils::MasterSlave::_slaveMode = false;
	utils::MasterSlave::_masterMode = false;
	utils::Parallel::clearGroups();
	utils::MasterSlave::_communication = nullptr;
}



}}} // namespace precice, geometry, tests

#endif // PRECICE_NO_MPI

