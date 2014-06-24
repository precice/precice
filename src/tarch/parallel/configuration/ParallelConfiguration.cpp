#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
#include "tarch/parallel/configuration/ParallelConfiguration.h"

#include "tarch/parallel/Node.h"


const std::string tarch::parallel::configuration::ParallelConfiguration::WaitWarningTime  = "time-out-warning";
const std::string tarch::parallel::configuration::ParallelConfiguration::DeadlockTimeout  = "deadlock-time-out";
const std::string tarch::parallel::configuration::ParallelConfiguration::Communicator     = "communicator";


tarch::logging::Log tarch::parallel::configuration::ParallelConfiguration::_log( "tarch::parallel::configuration::ParallelConfiguration" );


tarch::parallel::configuration::ParallelConfiguration::ParallelConfiguration():
  _timeoutWarning(-1),
  _deadlockTimeOut(-1),
#ifdef Parallel
  _communicator(MPI_COMM_WORLD),
#endif
  _validCommunicator(false) {
}


tarch::parallel::configuration::ParallelConfiguration::~ParallelConfiguration() {
}


bool tarch::parallel::configuration::ParallelConfiguration::writeWarning() const {
  return Node::getInstance().getRank()==0;
}


void tarch::parallel::configuration::ParallelConfiguration::parseSubtag( irr::io::IrrXMLReader* _xmlReader ) {
  assertion( _xmlReader != 0 );
  _validCommunicator = false;

  _timeoutWarning = _xmlReader->getAttributeValueAsInt( WaitWarningTime.c_str() );
  if ( _timeoutWarning==0 ) {
	if ( writeWarning() ) {
      _log.warning(
        "parseSubtag(...)",
        "switched off time-out warnings. Use attribute \"" + WaitWarningTime +
        "\" within tag <" + getTag() + "> to set time-out"
      );
    }
  }
  else if (_timeoutWarning*CLOCKS_PER_SEC < 1 ) {
	  if ( writeWarning() )
      _log.error(
        "parseSubgetTag()(...)", "attribute \"" + WaitWarningTime + "\" within getTag() <" + getTag() +
        "> has to have int value greater than 0 and must not create an overflow"
      );
    _timeoutWarning = -1;
  }

  _deadlockTimeOut = _xmlReader->getAttributeValueAsInt( DeadlockTimeout.c_str() );
  if ( _deadlockTimeOut==0 ) {
	if ( writeWarning() ) {
      _log.warning(
        "parseSubgetTag()(...)", "switched off deadlock time-out. Use attribute \"" + DeadlockTimeout +
        "\" within getTag() <" + getTag() + "> to set time-out"
      );
    }
  }
  else if (_deadlockTimeOut*CLOCKS_PER_SEC < 1 ) {
	  if ( writeWarning() )
      _log.error(
        "parseSubgetTag()(...)", "attribute \"" + DeadlockTimeout + "\" within getTag() <" + getTag() +
        "> has to have int value greater than 0 and must not create an overflow"
      );
    _deadlockTimeOut = -1;
  }


  if (_deadlockTimeOut<_timeoutWarning && _deadlockTimeOut>0) {
	if ( writeWarning() )
      _log.error(
        "parseSubgetTag()(...)", "attribute \"" + DeadlockTimeout +"\" within getTag() <" + getTag() +
        "> has to be greater/equal value of attribute timeout-warning"
      );
    _deadlockTimeOut = -1;
  }

  #ifdef Parallel
  if ( _xmlReader->getAttributeValue(Communicator.c_str())==0 ) {
    if ( writeWarning() ) {
        _log.error(
          "parseSubgetTag()(...)", "attribute \"" + Communicator +"\" within getTag() <" + getTag() +
          "> is missing"
        );
      _validCommunicator = false;
    }
  }
  else if ("default" == std::string(_xmlReader->getAttributeValue(Communicator.c_str()))) {
    _validCommunicator = true;
    _communicator      = MPI_COMM_WORLD;
  }
  else {
    if ( writeWarning() )
      _log.error(
        "parseSubgetTag()(...)", "attribute \"" + Communicator +"\" within getTag() <" + getTag() +
        "> has to have the value \"default\"; different communicators are not supported yet"
      );
    _validCommunicator = false;
  }
  #endif
}


bool tarch::parallel::configuration::ParallelConfiguration::isValid() const {
  bool valid = true;

  if(_timeoutWarning < 0) {
    logError("isValid", "attribute " << WaitWarningTime << " within tag <" + getTag() + "> has to be equal or larger than zero.");
    valid = false;
  }

  if(_deadlockTimeOut < 0) {
    logError("isValid", "attribute " << DeadlockTimeout << " within tag <" + getTag() + "> has to be equal or larger than zero.");
    valid = false;
  }

#ifdef Parallel
  if(!_validCommunicator) {
    logError("isValid", "attribute " << Communicator << " within tag <" + getTag() + "> has to be equal or larger than -1 but is " << _communicator << ".");
    valid = false;
  }
#endif

  return valid;
}


std::string tarch::parallel::configuration::ParallelConfiguration::getTag() const {
  return "parallel";
}


void tarch::parallel::configuration::ParallelConfiguration::toXML(std::ostream& out) const {
  out << "<!--" << std::endl
      << "  Configures the parallel node. " << std::endl
      << std::endl
      << "    | attribute name | semantics | allowed values " << std::endl
      << "    | time-out-warning  | After what time (ms) shall code write a warning. | Any possitive number of zero to switch off." << std::endl
      << "    | deadlock-time-out | After what time (ms) shall code stop with a time-out if no MPI message arrives. | Any possitive number of zero to switch off. " << std::endl
      << "    | communicator      | Which communicator shall the code use. | Any positive number or 'default'. " << std::endl
      << std::endl
      << "  -->" << std::endl;
  std::ostringstream result;
  result << "<" << getTag() << " "
         << WaitWarningTime             << "=\"" << _timeoutWarning << "\" "
         << DeadlockTimeout             << "=\"" << _deadlockTimeOut << "\" "
#ifdef Parallel
         << Communicator                << "=\"" << _communicator << "\" />"
#endif
         ;
}


void tarch::parallel::configuration::ParallelConfiguration::interpreteConfiguration() const{
  assertion(isValid());
  tarch::parallel::Node::getInstance().setTimeOutWarning(_timeoutWarning);
  tarch::parallel::Node::getInstance().setDeadlockTimeOut(_deadlockTimeOut);
#ifdef Parallel
  tarch::parallel::Node::getInstance().setCommunicator( _communicator );
#endif
}
