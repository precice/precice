#pragma once

#include "math/math.hpp"
#include "utils/Parallel.hpp"
#include <boost/test/unit_test.hpp>

#include "utils/MasterSlave.hpp"
#include "com/Communication.hpp"
#include "com/MPIDirectCommunication.hpp"
#include "m2n/M2N.hpp"
#include "m2n/GatherScatterComFactory.hpp"
#include "m2n/SharedPointer.hpp"


#ifndef PRECICE_NO_PETSC
#include "petscsys.h"
#endif

namespace precice
{
namespace testing
{

namespace bt = boost::unit_test;
using Par = precice::utils::Parallel;

/// Fixture to set and reset MPI communicator
struct MPICommRestrictFixture {
  explicit MPICommRestrictFixture(std::vector<int> &ranks)
  {
    // Restriction MUST always be called on all ranks, otherwise we hang
    if (static_cast<int>(ranks.size()) < Par::getCommunicatorSize()) {
      Par::restrictGlobalCommunicator(ranks);
    }
  }

  ~MPICommRestrictFixture()
  {
    Par::setGlobalCommunicator(Par::getCommunicatorWorld());
  }
};

/// Fixture to restrict to a single rank
/*
 * How does that differ from MPICommRestrictFixture({0})? The MPICommRestrictFixture restricts the communicator
 * to rank 0 and assigns that all ranks. This produces invalid communicators on all other ranks.
 * SingleRankFixture restricts every rank to itself. Effectively using MPI_COMM_SELF as communicator on each rank.
 * We don't use MPI_COMM_SELF because this causes errors when it's freed.
 */
struct SingleRankFixture {
  explicit SingleRankFixture()
  {
    // Restriction MUST always be called on all ranks, otherwise we hang
    Par::setGlobalCommunicator(Par::getRestrictedCommunicator({Par::getProcessRank()}));
  }

  ~SingleRankFixture()
  {
    Par::setGlobalCommunicator(Par::getCommunicatorWorld());
  }
};


/// Fixture to sync procceses before and after test
struct SyncProcessesFixture {
  SyncProcessesFixture()
  {
    Par::synchronizeProcesses();
  }

  ~SyncProcessesFixture()
  {
    Par::synchronizeProcesses();
  }
};

/// Fixture to create and destroy a master communication
/**
 * Many tests with parallel features need a working master communication. This fixture avoid code duplication for such tests.
 */
struct MasterComFixture {
  MasterComFixture()
  {
    utils::MasterSlave::_communication = com::PtrCommunication(new com::MPIDirectCommunication());
    int size = Par::getCommunicatorSize();

    if (utils::Parallel::getProcessRank() == 0){ //Master
      utils::Parallel::splitCommunicator( "Master" );
      utils::MasterSlave::_rank = 0;
      utils::MasterSlave::_size = size;
      utils::MasterSlave::_slaveMode = false;
      utils::MasterSlave::_masterMode = true;
      utils::MasterSlave::_communication->acceptConnection ( "Master", "Slaves", 0, 1);
      utils::MasterSlave::_communication->setRankOffset(1);
    }
    else {//Slaves
      assertion(utils::Parallel::getProcessRank() > 0 && utils::Parallel::getProcessRank() < size);
      utils::Parallel::splitCommunicator( "Slaves" );
      utils::MasterSlave::_rank = utils::Parallel::getProcessRank();
      utils::MasterSlave::_size = size;
      utils::MasterSlave::_slaveMode = true;
      utils::MasterSlave::_masterMode = false;
      utils::MasterSlave::_communication->requestConnection( "Master", "Slaves", utils::Parallel::getProcessRank()-1 , size-1 );
    }
  }

  ~MasterComFixture()
  {
    utils::MasterSlave::_communication = nullptr;
    utils::MasterSlave::reset();
    utils::Parallel::clearGroups();
  }
};

/// Fixture to create and destroy an m2n communication between two participants
struct M2NFixture {
  m2n::PtrM2N m2n;

  M2NFixture()
  {
    com::PtrCommunication participantCom =
          com::PtrCommunication(new com::MPIDirectCommunication());
    m2n::DistributedComFactory::SharedPointer distrFactory = m2n::DistributedComFactory::SharedPointer(
        new m2n::GatherScatterComFactory(participantCom));
    m2n = m2n::PtrM2N(new m2n::M2N(participantCom, distrFactory));

    if (utils::Parallel::getProcessRank() == 0){
      utils::Parallel::splitCommunicator( "ParticipantOne" );
      m2n->acceptMasterConnection ( "ParticipantOne", "ParticipantTwo");
    }
    else if(utils::Parallel::getProcessRank() == 1){//Master
      utils::Parallel::splitCommunicator( "ParticipantTwo" );
      m2n->requestMasterConnection ( "ParticipantOne", "ParticipantTwo" );
    }
  }

  ~M2NFixture()
  {
    utils::Parallel::clearGroups();
  }
};

/// Fixture to split two participants such that both can interact in an integration test
struct SplitParticipantsFixture {
  SplitParticipantsFixture()
  {
    if(utils::Parallel::getProcessRank()<=1){
      utils::Parallel::splitCommunicator( "ParticipantOne" );
      utils::Parallel::setGlobalCommunicator(utils::Parallel::getLocalCommunicator());
      utils::Parallel::clearGroups(); //This is important, if the testcase uses MPI communication again
    }
    else {
      assertion(utils::Parallel::getProcessRank() > 1 && utils::Parallel::getProcessRank() < 4);
      utils::Parallel::splitCommunicator( "ParticipantTwo" );
      utils::Parallel::setGlobalCommunicator(utils::Parallel::getLocalCommunicator());
      utils::Parallel::clearGroups();
    }
  }

  ~SplitParticipantsFixture()
  {
    utils::Parallel::setGlobalCommunicator(utils::Parallel::getCommunicatorWorld());
  }
};


/// Boost.Test decorator that makes the test run only on specfic ranks.
/*
 * This does not restrict the communicator, which must be done by installing the MPICommrestrictFixture.
 */
class OnRanks : public bt::decorator::base
{
public:
  explicit OnRanks(const std::vector<int> &ranks)
      : _ranks(ranks)
  {
  }

private:
  virtual void apply(bt::test_unit &tu)
  {
    size_t myRank = Par::getProcessRank();
    size_t size = Par::getCommunicatorSize();

    // If current rank is not in requested ranks
    if (std::find(_ranks.begin(), _ranks.end(), myRank) == _ranks.end()) {
      bt::framework::get<bt::test_suite>(tu.p_parent_id).remove(tu.p_id);
      return;
    }

    // If more ranks requested than available
    if (_ranks.size() > size) {
      bt::framework::get<bt::test_suite>(tu.p_parent_id).remove(tu.p_id);
      return;
    }

    // Install the fixture. Disabled because installing the fixture on just a
    // subset of ranks causes a restriction to be made from a subset of ranks
    // which means the application will hang.
    // tu.p_fixtures.value.push_back(
    //   bt::test_unit_fixture_ptr(
    //     new bt::class_based_fixture<MPICommRestrictFixture, std::vector<int>>(_ranks)));
  }

  virtual bt::decorator::base_ptr clone() const
  {
    return bt::decorator::base_ptr(new OnRanks(_ranks));
  }

  std::vector<int> _ranks;
};

/// Boost.Test decorator that makes the test run only on the master aka rank 0
class OnMaster : public OnRanks
{
public:
  explicit OnMaster()
      : OnRanks({0})
  {
  }
};

/// Boost.Test decorator that makes the test run only on a specific MPI size
class OnSize : public bt::decorator::base
{
public:
  explicit OnSize(const int size)
      : givenSize(size)
  {
  }

  virtual void apply(bt::test_unit &tu)
  {
    if (givenSize != Par::getCommunicatorSize()) {
      bt::framework::get<bt::test_suite>(tu.p_parent_id).remove(tu.p_id);
      return;
    }
  }

  virtual bt::decorator::base_ptr clone() const
  {
    return bt::decorator::base_ptr(new OnSize(givenSize));
  }

  const int givenSize;
};

/// Boost.Test decorator that deletes the test, unless a minimum number of ranks is available.
/*
 * This does not restrict the communicator, which must be done by installing the MPICommrestrictFixture.
 */
class MinRanks : public bt::decorator::base
{
public:
  explicit MinRanks(const int minimumSize)
      : minSize(minimumSize)
  {
  }

private:
  virtual void apply(bt::test_unit &tu)
  {
    if (minSize > Par::getCommunicatorSize()) {
      bt::framework::get<bt::test_suite>(tu.p_parent_id).remove(tu.p_id);
      return;
    }
  }

  virtual bt::decorator::base_ptr clone() const
  {
    return bt::decorator::base_ptr(new MinRanks(minSize));
  }

  const int minSize;
};


/// Boost.Test decorator that unconditionally deletes the test.
class Deleted : public bt::decorator::base
{
private:
  virtual void apply(bt::test_unit &tu)
  {
    bt::framework::get<bt::test_suite>(tu.p_parent_id).remove(tu.p_id);
  }

  virtual bt::decorator::base_ptr clone() const
  {
    return bt::decorator::base_ptr(new Deleted());
  }
};




/// equals to be used in tests. Prints both operatorans on failure
template <class DerivedA, class DerivedB>
boost::test_tools::predicate_result equals(const Eigen::MatrixBase<DerivedA> &A,
                                           const Eigen::MatrixBase<DerivedB> &B,
                                           double tolerance = math::NUMERICAL_ZERO_DIFFERENCE)
{
  if (not math::equals(A, B, tolerance)) {
    boost::test_tools::predicate_result res(false);
    Eigen::IOFormat format;
    if (A.cols() == 1) {
      format.rowSeparator = ", ";
      format.rowPrefix = "  ";
    }
    res.message() << "\n"
                  << A.format(format) << " != \n"
                  << B.format(format);
    return res;
  }
  return true;
}

/// equals to be used in tests. Prints both operatorans on failure
template <class Scalar>
typename std::enable_if<std::is_arithmetic<Scalar>::value, boost::test_tools::predicate_result>::type equals(const Scalar a, const Scalar b, const Scalar tolerance = math::NUMERICAL_ZERO_DIFFERENCE)
{
  if (not math::equals(a, b, tolerance)) {
    boost::test_tools::predicate_result res(false);
    res.message() << "Not equal: " << a << "!=" << b;
    return res;
  }
  return true;
}
}
} // namespace precice, testing
