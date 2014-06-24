#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
/*
 * ArgumentSet.cpp
 *
 *  Created on: 20.07.2011
 *      Author: michael
 */

#include "tarch/argument/ArgumentSet.h"
#include <iostream>
#include <sstream>
#include "tarch/Assertions.h"


tarch::logging::Log tarch::argument::ArgumentSet::_log("ArgumentSet()");

tarch::argument::ArgumentSet::ArgumentSet(std::string progId, int nrArgs):
_name(progId),
_nrArgs(nrArgs)
{
	// TODO Auto-generated constructor stub
}

tarch::argument::ArgumentSet::~ArgumentSet() {
	// TODO Auto-generated destructor stub
}

void tarch::argument::ArgumentSet::addArgument(std::string name, Argument::TYPE type, char* argument ) {
  assertion(_arguments.size()+1 <= _nrArgs);
	Argument newArgument(name, type, argument);
	_arguments.push_back( Argument(name, type, argument) );
}

void tarch::argument::ArgumentSet::addArgument(std::string name, Argument::TYPE type) {
 assertion(_arguments.size()+1 <= _nrArgs);
	_arguments.push_back( Argument(name, type) );
}

void tarch::argument::ArgumentSet::initialize(unsigned int argc, char* argv[]) {
	assertion2((argc-2)==_nrArgs, argc-2, _nrArgs);
	if((argc-2) != _nrArgs) {
		_log.error("initialize", "Error! Wrong number of arguments!");
		printDefaultArguments();
	}
	for(unsigned int i=0; i<_arguments.size(); i++) {
#ifdef Debug
		std::stringstream ss;
		ss << "\nArgument Nr. " <<i <<" : " <<argv[i+2];
		_log.debug("initialize",ss.str());
#endif
		_arguments[i].setArgument( std::string(argv[i+2]) );
	}
	_log.debug("initizalize()", "Argumentset initialization successful.");
}

double tarch::argument::ArgumentSet::getArgumentAsDouble( const std::string& argumentId ) const {
	return getArgument(argumentId).getArgumentAsDouble();
}

int    tarch::argument::ArgumentSet::getArgumentAsInt(const std::string& argumentId) const {
return getArgument(argumentId).getArgumentAsInt();
}

const char*  tarch::argument::ArgumentSet::getArgumentAsCharPointer(const std::string& argumentId) const {
return getArgument(argumentId).getArgumentAsCString();
}

std::string tarch::argument::ArgumentSet::getArgumentSetName() const {
	return _name;
}

bool tarch::argument::ArgumentSet::isArgumentSetName(const std::string& name) const {
	return (_name.compare(name) == 0);
}

tarch::argument::Argument tarch::argument::ArgumentSet::getArgument(const std::string argumentId) const {
	for (std::vector<Argument>::const_iterator it = _arguments.begin(); it!=_arguments.end(); ++it) {
		if((*it).isArgumentName(argumentId)) {
			return *it;
		}
	}
	std::stringstream ss;
	ss <<"\nERROR! Element " <<std::endl <<argumentId <<std::endl << "not found!!\n";
	print(ss);
	_log.error("getArgument", ss.str());
	exit(-1);
	return _arguments[0];
}

void tarch::argument::ArgumentSet::print() const {
	assertion(isValid());
	std::stringstream ss;
	ss <<"Program called:   " <<_name <<std::endl;
	for(unsigned int i=0; i<_arguments.size(); i++) {
	    ss <<_arguments[i].getName() <<": "<<_arguments[i].getArgumentAsString() <<std::endl;
	}
}

void tarch::argument::ArgumentSet::print(std::stringstream& ss) const {
	assertion(isValid());
	ss <<"Program called:   " <<_name <<std::endl;
	for(unsigned int i=0; i<_arguments.size(); i++) {
	    ss <<_arguments[i].getName() <<": " <<_arguments[i].getArgumentAsString() <<std::endl;
	}
}

void tarch::argument::ArgumentSet::printDefaultArguments() const {
	std::stringstream ss;
	print(ss);
	std::cout <<ss.str();
}

void tarch::argument::ArgumentSet::printDefaultArguments(std::stringstream& ss) const {
//	assertion(_arguments.size()>0);
	ss <<"Program:   " <<_name <<std::endl <<"   ";
	for(unsigned int i=0; i<_arguments.size(); i++) {
	    ss <<_arguments[i].getName() <<" ";
	}
}

bool tarch::argument::ArgumentSet::isValid() const {
	if(_nrArgs != _arguments.size()){
		_log.debug("isValid","_nrArgs invalid");
		return false;
	}

	for(unsigned int i=0; i<_arguments.size(); i++) {
		if(!_arguments[i].isValid()) {
			return false;
		}
	}
	return true;
}
