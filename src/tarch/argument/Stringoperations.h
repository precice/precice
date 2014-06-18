/*
 * Stringoperations.h
 *
 *  Created on: Mar 24, 2012
 *      Author: Michael Lieb
 *
 *  Some String operations which are used in different projects
 */

#ifndef STRINGOPERATIONS_H_
#define STRINGOPERATIONS_H_

namespace tarch {
  namespace stringoperations {
    /**
     * Replaces a the substring which begins with a string
     * given in variable from and ends with a string given
     * in the variable toFirstOccurenceOf with the string
     * given by replacement.
     *
     * Example usage:
     * string s("abcdefghijkl");
     * a = tarch::stringoperations::stringReplaceFromTo(
     *      s, "cd", "ij", "FOO");
     * a is now abFOOkl
     */
    inline std::string stringReplaceSubstringFromTo(
        std::string string,
        const std::string from,
        const std::string toFirstOccurenceOf,
        const std::string replacement) {
      return string.replace(
          string.find(from),
          (string.substr(string.find(from),string.size()-1)).find(toFirstOccurenceOf),
          replacement);
    }
  }
}




#endif /* STRINGOPERATIONS_H_ */
