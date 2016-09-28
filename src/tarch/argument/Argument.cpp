#ifndef PRECICE_NO_MPI
#include "mpi.h"
#endif
/*
 * Argument.cpp
 *
 *  Created on: 20.07.2011
 *      Author: michael
 */

#include "tarch/argument/Argument.h"
#include "utils/assertion.hpp"
#include <sstream>

logging::Logger tarch::argument::Argument::_log("Argument()");

tarch::argument::Argument::Argument(	std::string name, TYPE type, char* argument):
		_name(name),
		_type(type),
		_argument(NULL)
{
}

tarch::argument::Argument::Argument(	std::string name, TYPE type):
		_name(name),
		_type(type),
		_argument(const_cast<char*>(std::string("uninitialized").c_str()))
{
#ifdef Debug
	std::stringstream ss;
	ss << "Argument: " << name
	   << "Type    : " <<type;
	DEBUG("Argument",ss.str());
#endif
}

tarch::argument::Argument::~Argument() {
	// TODO Auto-generated destructor stub
}

double tarch::argument::Argument::string2double( const std::string& a ) const
{
	// Convert a string representation of a number into a floating point value.
	// Throws an int if the string contains anything but whitespace and a valid
	// numeric representation.
	//
	double result;
	std::string s( a );

	// Get rid of any trailing whitespace
	s.erase( s.find_last_not_of( " \f\n\r\t\v" ) + 1 );

	// Read it into the target type
	std::istringstream ss( s );
	ss >> result;

	// Check to see that there is nothing left over
	if (!ss.eof())
		throw 1;

	return result;
}

double tarch::argument::Argument::getArgumentAsDouble() const {
	double value = 0;
	try {
		value = string2double(_argument);
		//		 std:: cout <<"wunderbar double ist: " <<a;
	}
	catch (int i) {
		std::stringstream ss;
		ss << "Argument conversion in double failed: " << _argument;
		preciceWarnung("getArgumentAsDouble",ss.str());
	}
	return value;
}

bool tarch::argument::Argument::isArgInt(char const *c) const {
	while (*c) {
		if (!isdigit(*c++)) return false;
	}
	return true;
}

int tarch::argument::Argument::getArgumentAsInt() const {
	assertion(_argument!=NULL);
	assertion(isArgInt(_argument));
	return atoi(_argument);
}

const char* tarch::argument::Argument::getArgumentAsCString() const {
  const char* undefined = "undefined";
  return  (_argument == NULL) ? undefined :_argument;
}

std::string tarch::argument::Argument::getArgumentAsString() const {
	return (_argument == NULL) ? std::string("undefined") : std::string(_argument);
}

bool tarch::argument::Argument::isArgumentName(const std::string& name) const {
	return (_name == name);
}

std::string tarch::argument::Argument::getName() const {
	return _name;
}

void tarch::argument::Argument::setArgument(std::string argument) {
	char * buffer = new char[argument.length()];
	strcpy(buffer,argument.c_str());
	_argument = buffer;
#ifdef Debug
	std::stringstream ss;
	ss <<"Setting argument " << argument <<" successful\n";
	DEBUG("setArgument()", ss.str() );
#endif
}

bool tarch::argument::Argument::isValid() const {
	if(_name.empty()) {
		DEBUG("isValid", "_name invalid");
		return false;
	}
	if(_argument == NULL) {
		std::stringstream ss;
		ss << "In Argument " <<_name <<": _argument invalid";
		DEBUG("isValid", ss.str());
		return false;
	}
	if(_type != TYPE_INT && _type != TYPE_DOUBLE && _type != TYPE_STRING) {
		std::stringstream ss;
		ss << "In Argument " <<_name <<": _type invalid";
		DEBUG("isValid", ss.str());
		return false;
	}
	return true;
}
