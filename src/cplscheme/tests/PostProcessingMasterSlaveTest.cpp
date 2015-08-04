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
#include "cplscheme/ParallelCouplingScheme.hpp"
#include "cplscheme/impl/ConvergenceMeasure.hpp"
#include "cplscheme/impl/AbsoluteConvergenceMeasure.hpp"
#include "cplscheme/impl/MinIterationConvergenceMeasure.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "cplscheme/Constants.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "tarch/la/Vector.h"
#include "tarch/la/WrappedVector.h"

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
      testMethod (testVIQNPP);
      Par::setGlobalCommunicator(Par::getCommunicatorWorld());
    }
    //comm = Par::getRestrictedCommunicator(ranksWanted);
    //if (Par::getProcessRank() <= 4){
    //  Par::setGlobalCommunicator(comm); //necessary to be able to re-initialize with different leading ranks
    //  testMethod ( testMVQNPP() );
    //  Par::setGlobalCommunicator(Par::getCommunicatorWorld());
    //}
  }
}

void PostProcessingMasterSlaveTest::testVIQNPP()
{
	preciceTrace ( "testVIQNPP" ); assertion ( utils::Parallel::getCommunicatorSize() == 4 );

	//com::Communication::SharedPointer participantCom =
	//		com::Communication::SharedPointer(
	//				new com::MPIDirectCommunication());
	//m2n::DistributedComFactory::SharedPointer distrFactory =
	//		m2n::DistributedComFactory::SharedPointer(
	//				new m2n::GatherScatterComFactory(participantCom));
	//m2n::M2N::SharedPointer m2n = m2n::M2N::SharedPointer(
	//		new m2n::M2N(participantCom, distrFactory));
	com::Communication::SharedPointer masterSlaveCom =
			com::Communication::SharedPointer(
					new com::MPIDirectCommunication());
	utils::MasterSlave::_communication = masterSlaveCom;

	utils::Parallel::synchronizeProcesses();

	//if (utils::Parallel::getProcessRank() == 0) { //NASTIN
	//	utils::Parallel::splitCommunicator("NASTIN");
	//	m2n->acceptMasterConnection("NASTIN", "SOLIDZMaster");
	//} else
	if (utils::Parallel::getProcessRank() == 0) { //Master
		utils::Parallel::splitCommunicator("SOLIDZMaster");
		//m2n->requestMasterConnection("NASTIN", "SOLIDZMaster");
	} else if (utils::Parallel::getProcessRank() == 1) { //Slave1
		utils::Parallel::splitCommunicator("SOLIDZSlaves");
	} else if (utils::Parallel::getProcessRank() == 2) { //Slave2
		utils::Parallel::splitCommunicator("SOLIDZSlaves");
	} else if (utils::Parallel::getProcessRank() == 3) { //Slave3
		utils::Parallel::splitCommunicator("SOLIDZSlaves");
	}

	if (utils::Parallel::getProcessRank() == 0) { //Master
		masterSlaveCom->acceptConnection("SOLIDZMaster", "SOLIDZSlaves", 0, 1);
		masterSlaveCom->setRankOffset(0);
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
	double singularityLimit = 1e-10;
	std::vector<int> dataIDs;
	dataIDs.push_back(0);
	dataIDs.push_back(1);
	std::map<int, double> scalings;
	scalings.insert(std::make_pair(0,1.0));
	scalings.insert(std::make_pair(1,1.0));
	mesh::PtrMesh dummyMesh ( new mesh::Mesh("dummyMesh", 3, false) );

	cplscheme::impl::IQNILSPostProcessing pp(initialRelaxation,maxIterationsUsed,
										   timestepsReused, singularityLimit, dataIDs, scalings);

	utils::DynVector dvalues;
	utils::DynVector dcol1;
	utils::DynVector fvalues;
	utils::DynVector fcol1;

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
		dvalues.append(1.0); dvalues.append(2.0); dvalues.append(3.0); dvalues.append(4.0);
		dcol1.append(1.0); dcol1.append(1.0); dcol1.append(1.0); dcol1.append(1.0);

		PtrCouplingData dpcd(new CouplingData(&dvalues,dummyMesh,false,1));

		//init forces
		fvalues.append(0.1); fvalues.append(0.1); fvalues.append(0.1); fvalues.append(0.1);
		fcol1.append(0.2); fcol1.append(0.2); fcol1.append(0.2); fcol1.append(0.2);

		PtrCouplingData fpcd(new CouplingData(&fvalues,dummyMesh,false,1));

		data.insert(std::pair<int,PtrCouplingData>(0,dpcd));
		data.insert(std::pair<int,PtrCouplingData>(1,fpcd));

		pp.initialize(data);

		dpcd->oldValues.column(0) = dcol1;
		fpcd->oldValues.column(0) = fcol1;

	} else if (utils::Parallel::getProcessRank() == 1) { //Slave1
		utils::MasterSlave::_rank = 1;
		utils::MasterSlave::_size = 4;
		utils::MasterSlave::_slaveMode = true;
		utils::MasterSlave::_masterMode = false;

		/**
		 * processor with 4 vertices
		 */

		//init displacements
		dvalues.append(5.0); dvalues.append(6.0); dvalues.append(7.0); dvalues.append(8.0);
		dcol1.append(1.0); dcol1.append(1.0); dcol1.append(1.0); dcol1.append(1.0);

		PtrCouplingData dpcd(new CouplingData(&dvalues,dummyMesh,false,1));

		//init forces
		fvalues.append(0.1); fvalues.append(0.1); fvalues.append(0.1); fvalues.append(0.1);
		fcol1.append(0.2); fcol1.append(0.2); fcol1.append(0.2); fcol1.append(0.2);

		PtrCouplingData fpcd(new CouplingData(&fvalues,dummyMesh,false,1));

		data.insert(std::pair<int,PtrCouplingData>(0,dpcd));
		data.insert(std::pair<int,PtrCouplingData>(1,fpcd));

		pp.initialize(data);

		dpcd->oldValues.column(0) = dcol1;
		fpcd->oldValues.column(0) = fcol1;

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

		dpcd->oldValues.column(0) = dcol1;
		fpcd->oldValues.column(0) = fcol1;

	} else if (utils::Parallel::getProcessRank() == 3) { //Slave3
		utils::MasterSlave::_rank = 3;
		utils::MasterSlave::_size = 4;
		utils::MasterSlave::_slaveMode = true;
		utils::MasterSlave::_masterMode = false;

		/**
		 * processor with 2 vertices
		 */

		//init displacements
		dvalues.append(1.0); dvalues.append(2.0);
		dcol1.append(1.0); dcol1.append(1.0);

		PtrCouplingData dpcd(new CouplingData(&dvalues,dummyMesh,false,1));

		//init forces
		fvalues.append(0.1); fvalues.append(0.1);
		fcol1.append(0.2); fcol1.append(0.2);

		PtrCouplingData fpcd(new CouplingData(&fvalues,dummyMesh,false,1));

		data.insert(std::pair<int,PtrCouplingData>(0,dpcd));
		data.insert(std::pair<int,PtrCouplingData>(1,fpcd));

		pp.initialize(data);

		dpcd->oldValues.column(0) = dcol1;
		fpcd->oldValues.column(0) = fcol1;
	}

	std::cout<<"perform pp"<<std::endl;


	pp.performPostProcessing(data);

	std::cout<<"performed"<<std::endl;



	std::cout<<"First Iteration:"<<std::endl;
	/*
	validateWithParams1(tarch::la::equals((*data.at(0)->values)(0), 1.00), (*data.at(0)->values)(0));
	validateWithParams1(tarch::la::equals((*data.at(0)->values)(1), 1.01), (*data.at(0)->values)(1));
	validateWithParams1(tarch::la::equals((*data.at(0)->values)(2), 1.02), (*data.at(0)->values)(2));
	validateWithParams1(tarch::la::equals((*data.at(0)->values)(3), 1.03), (*data.at(0)->values)(3));
	validateWithParams1(tarch::la::equals((*data.at(1)->values)(0), 0.199), (*data.at(1)->values)(0));
	validateWithParams1(tarch::la::equals((*data.at(1)->values)(1), 0.199), (*data.at(1)->values)(1));
	validateWithParams1(tarch::la::equals((*data.at(1)->values)(2), 0.199), (*data.at(1)->values)(2));
	*/

	utils::DynVector newdvalues;
	if (utils::Parallel::getProcessRank() == 0) { //Master

		std::cout<<"  Master:"<<std::endl;
		std::cout<<"(*data.at(0)->values)(0): "<<(*data.at(0)->values)(0)<<std::endl;
		std::cout<<"(*data.at(0)->values)(1): "<<(*data.at(0)->values)(1)<<std::endl;
		std::cout<<"(*data.at(0)->values)(2): "<<(*data.at(0)->values)(2)<<std::endl;
		std::cout<<"(*data.at(0)->values)(3): "<<(*data.at(0)->values)(3)<<std::endl;
		std::cout<<"(*data.at(1)->values)(0): "<<(*data.at(1)->values)(0)<<std::endl;
		std::cout<<"(*data.at(1)->values)(1): "<<(*data.at(1)->values)(1)<<std::endl;
		std::cout<<"(*data.at(1)->values)(2): "<<(*data.at(1)->values)(2)<<std::endl;
		std::cout<<"(*data.at(1)->values)(3): "<<(*data.at(1)->values)(3)<<std::endl;

		newdvalues.append(10.0); newdvalues.append(10.0); newdvalues.append(10.0); newdvalues.append(10.0);

	} else if (utils::Parallel::getProcessRank() == 1) { //Slave1

		std::cout<<"  Slave 1:"<<std::endl;
		std::cout<<"(*data.at(0)->values)(0): "<<(*data.at(0)->values)(0)<<std::endl;
		std::cout<<"(*data.at(0)->values)(1): "<<(*data.at(0)->values)(1)<<std::endl;
		std::cout<<"(*data.at(0)->values)(2): "<<(*data.at(0)->values)(2)<<std::endl;
		std::cout<<"(*data.at(0)->values)(3): "<<(*data.at(0)->values)(3)<<std::endl;
		std::cout<<"(*data.at(1)->values)(0): "<<(*data.at(1)->values)(0)<<std::endl;
		std::cout<<"(*data.at(1)->values)(1): "<<(*data.at(1)->values)(1)<<std::endl;
		std::cout<<"(*data.at(1)->values)(2): "<<(*data.at(1)->values)(2)<<std::endl;
		std::cout<<"(*data.at(1)->values)(3): "<<(*data.at(1)->values)(3)<<std::endl;

		newdvalues.append(10.0); newdvalues.append(10.0); newdvalues.append(10.0); newdvalues.append(10.0);

	} else if (utils::Parallel::getProcessRank() == 2) { //Slave2
		// empty proc

	} else if (utils::Parallel::getProcessRank() == 3) { //Slave3

		std::cout<<"  Slave 3:"<<std::endl;
		std::cout<<"(*data.at(0)->values)(0): "<<(*data.at(0)->values)(0)<<std::endl;
		std::cout<<"(*data.at(0)->values)(1): "<<(*data.at(0)->values)(1)<<std::endl;
		std::cout<<"(*data.at(1)->values)(0): "<<(*data.at(1)->values)(0)<<std::endl;
		std::cout<<"(*data.at(1)->values)(1): "<<(*data.at(1)->values)(1)<<std::endl;

		newdvalues.append(10.0); newdvalues.append(10.0);
	}

	data.begin()->second->values = &newdvalues;

	pp.performPostProcessing(data);

	/*
	validateWithParams1(tarch::la::equals((*data.at(0)->values)(0), -5.63855295490201413600e-01), (*data.at(0)->values)(0));
	validateWithParams1(tarch::la::equals((*data.at(0)->values)(1), 6.09906404008709657205e-01), (*data.at(0)->values)(1));
	validateWithParams1(tarch::la::equals((*data.at(0)->values)(2), 1.78366810350762072801e+0), (*data.at(0)->values)(2));
	validateWithParams1(tarch::la::equals((*data.at(0)->values)(3), 2.95742980300653179881e+00), (*data.at(0)->values)(3));
	validateWithParams1(tarch::la::equals((*data.at(1)->values)(0), 8.27975917496077823410e-02), (*data.at(1)->values)(0));
	validateWithParams1(tarch::la::equals((*data.at(1)->values)(1), 8.27975917496077823410e-02), (*data.at(1)->values)(1));
	validateWithParams1(tarch::la::equals((*data.at(1)->values)(2), 8.27975917496077823410e-02), (*data.at(1)->values)(2));
	*/


	utils::MasterSlave::_slaveMode = false;
	utils::MasterSlave::_masterMode = false;
	utils::Parallel::synchronizeProcesses();
	utils::Parallel::clearGroups();
}



}}} // namespace precice, geometry, tests

#endif // PRECICE_NO_MPI

