/*
 * ArgumentSetFabric.h
 *
 *  Created on: 21.07.2011
 *      Author: michael
 */

#ifndef ArgumentSetFabric_H_
#define ArgumentSetFabric_H_

#include "tarch/argument/ArgumentSet.h"
#include <vector>
#include "tarch/logging/Log.h"
namespace tarch {
  namespace argument {
  class ArgumentSetFabric;
}}

class tarch::argument::ArgumentSetFabric {
public:
  ArgumentSetFabric();
  static ArgumentSet getArgumentSet(unsigned int argc, char* argv[]);
  static void print();
  static void print(std::stringstream& ss);
  static void printDefaultArguments(std::stringstream& ss);
  static void printDefaultArguments();
  static void addArgumentSet(ArgumentSet argumentSet);
  static bool isCurrentSet(std::string set);
  virtual ~ArgumentSetFabric();
  static void printAllPossibleConfigurationsAndExit();
//  static std::pair<std::string,unsigned int> getCurrentSetting(unsigned int argc, char* argv[]);

private:
  static ArgumentSet getSpecificSet(unsigned int argc, char* argv[]);

private:
  static tarch::logging::Log _log;
  static std::vector<tarch::argument::ArgumentSet> _argumentSets;
  static bool _isInitialized;
};

#endif /* ArgumentSetFabric_H_ */
