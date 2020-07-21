#pragma once

#include <algorithm>
#include <exception>
#include <initializer_list>
#include <map>
#include <sstream>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>
#include "logging/Logger.hpp"
#include "utils/assertion.hpp"
#include "xml/ValueParser.hpp"

namespace precice {
namespace xml {

template <typename ATTRIBUTE_T>
class XMLAttribute {
  static_assert(std::is_default_constructible<ATTRIBUTE_T>::value, "The value type of XMLAttributes need to be default-constructible.");

public:
  XMLAttribute() = delete;

  explicit XMLAttribute(std::string name)
      : _name(std::move(name)){};

  XMLAttribute(std::string name, ATTRIBUTE_T defaultValue)
      : _name(std::move(name)), _hasDefaultValue(true), _defaultValue(std::move(defaultValue)){};

  XMLAttribute(const XMLAttribute<ATTRIBUTE_T> &other) = default;

  XMLAttribute &operator=(const XMLAttribute<ATTRIBUTE_T> &other) = default;

  /// Sets a documentation string for the attribute.
  XMLAttribute &setDocumentation(std::string documentation);

  const std::string &getUserDocumentation() const
  {
    return _doc;
  }

  XMLAttribute &setOptions(std::vector<ATTRIBUTE_T> options);

  template <class T>
  XMLAttribute &setOptions(std::initializer_list<T> &&options)
  {
    static_assert(std::is_convertible<T, ATTRIBUTE_T>::value, "Type of initializer_list must be converible to ATTRIBUTE_T!");
    return setOptions(std::vector<ATTRIBUTE_T>(options.begin(), options.end()));
  }

  const std::vector<ATTRIBUTE_T> &getOptions() const
  {
    return _options;
  };

  XMLAttribute &setDefaultValue(const ATTRIBUTE_T &defaultValue);

  const ATTRIBUTE_T &getDefaultValue() const
  {
    return _defaultValue;
  };

  bool hasDefaultValue() const
  {
    return _hasDefaultValue;
  };

  bool hasValidation() const
  {
    return _hasValidation;
  };

  void readValue(const std::map<std::string, std::string> &aAttributes);

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

private:
  logging::Logger _log{"xml::XMLAttribute"};

  std::string _name;

  std::string _doc;

  bool _read = false;

  ATTRIBUTE_T _value{};

  bool _hasDefaultValue = false;

  ATTRIBUTE_T _defaultValue{};

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
XMLAttribute<ATTRIBUTE_T> &XMLAttribute<ATTRIBUTE_T>::setDocumentation(std::string documentation)
{
  _doc = std::move(documentation);
  return *this;
}

template <typename ATTRIBUTE_T>
XMLAttribute<ATTRIBUTE_T> &XMLAttribute<ATTRIBUTE_T>::setOptions(std::vector<ATTRIBUTE_T> options)
{
  const auto iter = std::unique(options.begin(), options.end());
  _options        = std::vector<ATTRIBUTE_T>(options.begin(), iter);
  _hasValidation  = true;
  return *this;
}

template <typename ATTRIBUTE_T>
XMLAttribute<ATTRIBUTE_T> &XMLAttribute<ATTRIBUTE_T>::setDefaultValue(const ATTRIBUTE_T &defaultValue)
{
  PRECICE_TRACE(defaultValue);
  _hasDefaultValue = true;
  set(_defaultValue, defaultValue);
  return *this;
}

template <typename ATTRIBUTE_T>
void XMLAttribute<ATTRIBUTE_T>::readValue(const std::map<std::string, std::string> &aAttributes)
{
  PRECICE_TRACE(_name);
  PRECICE_ASSERT(!_read, "Attribute \"" + _name + "\" has already been read.");

  const auto position = aAttributes.find(getName());
  if (position == aAttributes.end()) {
    if (not _hasDefaultValue) {
      PRECICE_ERROR("Attribute \"" + _name + "\" is required, but was not defined.");
    }
    set(_value, _defaultValue);
  } else {
    try {
      readValueSpecific(position->second, _value);
    } catch (const std::exception &e) {
      PRECICE_ERROR(e.what());
    }
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
        for (; first != _options.end(); ++first) {
          stream << " or value must be \"" << *first << '"';
        }

        PRECICE_ERROR(stream.str());
      }
    }
  }
  PRECICE_DEBUG("Read valid attribute \"" << getName() << "\" value = " << _value);
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
inline XMLAttribute<std::string> makeXMLAttribute(std::string name, const char *defaultValue)
{
  return XMLAttribute<std::string>(std::move(name), defaultValue);
}

/** creates an XMLAttribute given a name and a default value.
 *  
 *  @param[in] name the name of the attribute
 *  @param[in] defaultValue the default value of the attribute
 *  @return an XMLAttribute with the above settings
 */
template <typename T>
XMLAttribute<T> makeXMLAttribute(std::string name, T defaultValue)
{
  return XMLAttribute<T>(std::move(name), std::move(defaultValue));
}

} // namespace xml
} // namespace precice
