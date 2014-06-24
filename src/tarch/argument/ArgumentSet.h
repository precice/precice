/*
 * ArgumentSet.h
 *
 *  Created on: 20.07.2011
 *      Author: michael
 */

#ifndef ARGUMENTSET_H_
#define ARGUMENTSET_H_

#ifdef Parallel
#include <mpi.h>
#endif
#include <vector>
#include <string>
#include <sstream>
#include "tarch/argument/Argument.h"
#include "tarch/logging/Log.h"

namespace tarch {
 namespace argument {
  class ArgumentSet;
}}

class tarch::argument::ArgumentSet {
public:
	ArgumentSet(std::string type, int nrArgs);
	virtual ~ArgumentSet();

	void initialize(unsigned int argc, char* argv[]);

	double getArgumentAsDouble(     const std::string& argumentId) const;
	int    getArgumentAsInt(        const std::string& argumentId) const;
	const char*  getArgumentAsCharPointer(const std::string& argumentId) const;

	std::string getArgumentSetName() const;
	bool isArgumentSetName(const std::string& name) const;

	void addArgument(std::string name, Argument::TYPE, char* argument);
	void addArgument(std::string name, Argument::TYPE);

	void print() const;
	void print(std::stringstream& ss) const;
	void printDefaultArguments() const;
	void printDefaultArguments(std::stringstream& ss) const;
	bool isValid() const;

private:
	static tarch::logging::Log _log;
	tarch::argument::Argument getArgument(std::string argumentId) const;
	std::string             _name;
	unsigned int            _nrArgs;
	std::vector< Argument > _arguments;
};

#endif /* ARGUMENTSET_H_ */
