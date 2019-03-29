#pragma once

#include <iostream>
#include <string>
#include <vector>
#include <type_traits>
#include <initializer_list>

#include "logging/Logger.hpp"
#include "math/math.hpp"
#include "utils/TypeNames.hpp"
#include "utils/assertion.hpp"
#include "utils/String.hpp"

namespace precice
{
namespace xml
{

template <typename ATTRIBUTE_T>
class XMLAttribute
{
public:
  XMLAttribute() = delete;

  explicit XMLAttribute(std::string name) : _name(std::move(name)) {};

  XMLAttribute(std::string name, ATTRIBUTE_T defaultValue): _name(std::move(name)), _hasDefaultValue(true), _defaultValue(std::move(defaultValue)) {};

  XMLAttribute(const XMLAttribute<ATTRIBUTE_T> &other) = default;

  XMLAttribute& operator=(const XMLAttribute<ATTRIBUTE_T> &other) = default;

  /// Sets a documentation string for the attribute.
  XMLAttribute& setDocumentation(std::string documentation);

  const std::string &getUserDocumentation() const
  {
    return _doc;
  }

  XMLAttribute& setOptions(std::vector<ATTRIBUTE_T> options);


  template<class T>
  XMLAttribute& setOptions(std::initializer_list<T>&& options)
  {
    static_assert(std::is_convertible<T, ATTRIBUTE_T>::value, "Type of initializer_list must be converible to ATTRIBUTE_T!");
    return setOptions(std::vector<ATTRIBUTE_T>(options.begin(), options.end()));
  }

  XMLAttribute& setDefaultValue(const ATTRIBUTE_T &defaultValue);

  void readValue(std::map<std::string, std::string> &aAttributes);

  void readValueSpecific(std::string &rawValue, double &value);

  void readValueSpecific(std::string &rawValue, int &value);

  void readValueSpecific(std::string &rawValue, std::string &value);

  void readValueSpecific(std::string &rawValue, bool &value);

  void readValueSpecific(std::string &rawValue, Eigen::VectorXd &value);

  //Eigen::VectorXd getAttributeValueAsEigenVectorXd(std::string& rawValue);

  const std::string &getName() const
  {
    return _name;
  };

  const ATTRIBUTE_T &getValue() const
  {
    return _value;
  };

  void setRead(bool read)
  {
    _read = read;
  };

  bool isRead() const
  {
    return _read;
  };

  /// Returns a documentation string about the attribute.
  std::string printDocumentation() const;

  std::string printDTD(const std::string &ElementName) const;

private:
  logging::Logger _log{"xml::XMLAttribute"};

  std::string _name;

  std::string _doc;

  bool _read = false;

  ATTRIBUTE_T _value;

  bool _hasDefaultValue = false;

  ATTRIBUTE_T _defaultValue;

  bool _hasValidation = false;

  std::vector<ATTRIBUTE_T> _options;

  /// Sets non Eigen::VectorXd type values.
  template <typename VALUE_T>
  typename std::enable_if<
      std::is_same<VALUE_T, ATTRIBUTE_T>::value && not std::is_same<VALUE_T, Eigen::VectorXd>::value, void>::type
  set(ATTRIBUTE_T &toSet, const VALUE_T &setter);

  /// Sets Eigen::VectorXd type values by clearing and copy.
  template <typename VALUE_T>
  typename std::enable_if<
      std::is_same<VALUE_T, ATTRIBUTE_T>::value && std::is_same<VALUE_T, Eigen::VectorXd>::value, void>::type
  set(ATTRIBUTE_T &toSet, const VALUE_T &setter);
};

template <typename ATTRIBUTE_T>
XMLAttribute<ATTRIBUTE_T>& XMLAttribute<ATTRIBUTE_T>::setDocumentation(std::string documentation)
{
  _doc = std::move(documentation);
  return *this;
}

template <typename ATTRIBUTE_T>
XMLAttribute<ATTRIBUTE_T>& XMLAttribute<ATTRIBUTE_T>::setOptions(std::vector<ATTRIBUTE_T> options)
{
  const auto iter = std::unique(options.begin(), options.end());
  _options     = std::vector<ATTRIBUTE_T>(options.begin(), iter);
  _hasValidation = true;
  return *this;
}

template <typename ATTRIBUTE_T>
XMLAttribute<ATTRIBUTE_T>& XMLAttribute<ATTRIBUTE_T>::setDefaultValue(const ATTRIBUTE_T &defaultValue)
{
  TRACE(defaultValue);
  _hasDefaultValue = true;
  set(_defaultValue, defaultValue);
  return *this;
}

template <typename ATTRIBUTE_T>
void XMLAttribute<ATTRIBUTE_T>::readValue(std::map<std::string, std::string> &aAttributes)
{
  TRACE(_name);
  if (_read) {
    std::cout << "Attribute \"" + _name + "\" is defined multiple times\n";
    ERROR("Attribute \"" + _name + "\" is defined multiple times");
  }

  if (aAttributes.find(getName()) == aAttributes.end()) {
    if (not _hasDefaultValue) {
      std::cout << "Attribute \"" + _name + "\" missing\n";
      ERROR("Attribute \"" + _name + "\" missing");
    }
    set(_value, _defaultValue);
  } else {
    readValueSpecific(aAttributes[getName()], _value);
    if (_hasValidation) {
      if (std::find(_options.begin(), _options.end(), _value) == _options.end()) {
        std::ostringstream stream;
        stream << "Invalid value \"" << _value << "\" of attribute \""
               << getName() << "\": ";
        // print first
        auto first = _options.begin();
        stream << "value must be \"" << *first << '"';
        ++first;
        // print the remaining with separator
        for(;first != _options.end();++first) {
            stream << " or value must be \"" << *first << '"';
        }

        std::cout << stream.str() << '\n';
        ERROR(stream.str());
      }
    }
  }
  DEBUG("Read valid attribute \"" << getName() << "\" value = " << _value);
}

template <typename ATTRIBUTE_T>
void XMLAttribute<ATTRIBUTE_T>::readValueSpecific(std::string &rawValue, double &value)
{
  try {
    if (rawValue.find('/') != std::string::npos) {
      std::string left  = rawValue.substr(0, rawValue.find('/'));
      std::string right = rawValue.substr(rawValue.find('/') + 1, rawValue.size() - rawValue.find('/') - 1);

      value = std::stod(left) / std::stod(right);
    } else {
      value = std::stod(rawValue);
    }
  } catch (...) {
    ERROR("String to Double error");
  }
}

template <typename ATTRIBUTE_T>
void XMLAttribute<ATTRIBUTE_T>::readValueSpecific(std::string &rawValue, int &value)
{
  try {
    value = std::stoi(rawValue);
  } catch (...) {
    ERROR("String to Int error");
  }
}

template <typename ATTRIBUTE_T>
void XMLAttribute<ATTRIBUTE_T>::readValueSpecific(std::string &rawValue, std::string &value)
{
  value = rawValue;
}

template <typename ATTRIBUTE_T>
void XMLAttribute<ATTRIBUTE_T>::readValueSpecific(std::string &rawValue, bool &value)
{
  value = precice::utils::convertStringToBool(rawValue);
}

template <typename ATTRIBUTE_T>
void XMLAttribute<ATTRIBUTE_T>::readValueSpecific(std::string &rawValue, Eigen::VectorXd &value)
{
  Eigen::VectorXd vec;

  std::string valueString(rawValue);
  bool        componentsLeft = true;
  int         i              = 0;
  while (componentsLeft) {
    std::string tmp1(rawValue);
    // erase entries before i-th entry
    for (int j = 0; j < i; j++) {
      if (tmp1.find(';') != std::string::npos) {
        tmp1.erase(0, tmp1.find(';') + 1);
      } else {
        componentsLeft = false;
      }
    }
    // if we are not in the last vector component...
    if (tmp1.find(';') != std::string::npos) {
      // ..., erase entries after i-th entry
      tmp1.erase(tmp1.find(';'), tmp1.size());
    }

    if (componentsLeft) {

      vec.conservativeResize(vec.rows() + 1);
      vec(vec.rows() - 1) = std::stod(tmp1);
    }
    i++;
  }

  value = vec;
}

/*template<>
Eigen::VectorXd XMLAttribute<Eigen::VectorXd>::getAttributeValueAsEigenVectorXd(std::string& rawValue)
{
  Eigen::VectorXd vec;

  std::string valueString(rawValue);
  bool componentsLeft = true;
  int i = 0;
  while (componentsLeft){
	std::string tmp1(rawValue);
	// erase entries before i-th entry
	for (int j = 0; j < i; j++){
	  if (tmp1.find(';') != std::string::npos){
		tmp1.erase(0,tmp1.find(';')+1);
	  }
	  else {
		componentsLeft = false;
	  }
	}
	// if we are not in the last vector component...
	if (tmp1.find(';') != std::string::npos){
	  // ..., erase entries after i-th entry
	  tmp1.erase(tmp1.find(';'),tmp1.size());
	}
	if (componentsLeft){
	   
	  vec.conservativeResize(vec.rows()+1);
	  vec(vec.rows()-1) = std::stod(tmp1);
	}
	i++;
  }
  return vec;
}*/

template <typename ATTRIBUTE_T>
std::string XMLAttribute<ATTRIBUTE_T>::printDTD(const std::string &ElementName) const
{
  std::ostringstream dtd;
  dtd << "<!ATTLIST " << ElementName << " " << _name << " CDATA ";

  if (_hasDefaultValue) {
    dtd << "\"" << _defaultValue << "\"";
  } else {
    dtd << "#REQUIRED";
  }

  dtd << ">\n";

  return dtd.str();
}

template <typename ATTRIBUTE_T>
std::string XMLAttribute<ATTRIBUTE_T>::printDocumentation() const
{
  std::ostringstream doc;
  doc << _name << "=\"{" << utils::getTypeName(_value);
  if (_hasValidation) {
    assertion(!_options.empty());
    doc << ":";
    // print the first item
    auto first = _options.begin();
    doc << '\'' << *first << '\'';
    ++first;
    // print the remaining items with separator
    for(;first != _options.end(); ++first) {
        doc << " or '" << *first << '\'';
    }
  }
  doc << "}";
  if (_hasDefaultValue) {
    doc << "(default:'" << _defaultValue << "')";
  }
  doc << "\"";
  return doc.str();
}

template <typename ATTRIBUTE_T>
template <typename VALUE_T>
typename std::enable_if<
    std::is_same<VALUE_T, ATTRIBUTE_T>::value && not std::is_same<VALUE_T, Eigen::VectorXd>::value, void>::type
XMLAttribute<ATTRIBUTE_T>::set(
    ATTRIBUTE_T &  toSet,
    const VALUE_T &setter)
{
  toSet = setter;
}

template <typename ATTRIBUTE_T>
template <typename VALUE_T>
typename std::enable_if<
    std::is_same<VALUE_T, ATTRIBUTE_T>::value && std::is_same<VALUE_T, Eigen::VectorXd>::value, void>::type
XMLAttribute<ATTRIBUTE_T>::set(
    ATTRIBUTE_T &  toSet,
    const VALUE_T &setter)
{
  toSet = setter;
}

/** creates an XMLAttribute given a name and a default value.
 *  
 *  @param[in] name the name of the attribute
 *  @param[in] defaultValue the default value of the attribute
 *  @return an XMLAttribute with the above settings
 */
inline XMLAttribute<std::string> makeXMLAttribute(std::string name, const char * defaultValue) {
    return XMLAttribute<std::string>(std::move(name), defaultValue);
}

/** creates an XMLAttribute given a name and a default value.
 *  
 *  @param[in] name the name of the attribute
 *  @param[in] defaultValue the default value of the attribute
 *  @return an XMLAttribute with the above settings
 */
template<typename T>
XMLAttribute<T> makeXMLAttribute(std::string name, T defaultValue) {
    return XMLAttribute<T>(std::move(name), std::move(defaultValue));
}

}
} // namespace precice, xml
