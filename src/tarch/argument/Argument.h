/*
 * Argument.h
 *
 *  Created on: 20.07.2011
 *      Author: michael
 */

#ifndef ARGUMENT_H_
#define ARGUMENT_H_

#ifdef Parallel
#include <mpi.h>
#endif
#include <cstdlib>
#include <cstring>
#include <string>
#include "tarch/logging/Log.h"

namespace tarch {
  namespace argument {
  class Argument;
}}

class tarch::argument::Argument {
public:
	enum TYPE {TYPE_DOUBLE, TYPE_INT, TYPE_STRING};

public:
	Argument( std::string name, TYPE type, char* argument );
	Argument( std::string name, TYPE type );
	virtual ~Argument();

	double getArgumentAsDouble() const;
	int getArgumentAsInt() const;
	const char* getArgumentAsCString() const;
	std::string getArgumentAsString() const;
	bool isArgumentName(const std::string& name) const;
	std::string getName() const;
	void setArgument(std::string argument);
	bool isArgInt(char const *arg) const;
	double string2double( const std::string& a ) const;
	bool isValid() const;

private:
	static tarch::logging::Log _log;
	std::string _name;
	TYPE        _type;
	char*       _argument;
};

#endif /* ARGUMENT_H_ */
