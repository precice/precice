// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_UTILS_XMLATTRIBUTE_HPP_
#define PRECICE_UTILS_XMLATTRIBUTE_HPP_

#include <type_traits>

#include "Validator.hpp"
#include "tarch/irr/XML.h"
#include "tarch/logging/Log.h"
#include "utils/Dimensions.hpp"
#include "utils/Helpers.hpp"
#include "utils/Globals.hpp"

#include <string>

namespace precice {
namespace utils {



template<typename ATTRIBUTE_T>
class XMLAttribute
{
public:

  typedef tarch::irr::io::IrrXMLReader XMLReader;

  /**
   * @brief Constructor for compatibility with map::operator[]. Not to be used!
   *
   * Gives an assertion on use.
   */
  XMLAttribute();

  XMLAttribute ( const std::string& name );

  XMLAttribute ( const XMLAttribute<ATTRIBUTE_T>& rhs );

  virtual ~XMLAttribute() { delete _validator; };

  /**
   * @brief Sets a documentation string for the attribute.
   */
  void setDocumentation ( const std::string& documentation );

  const std::string& getUserDocumentation() const
  {
    return _doc;
  }

  void setValidator ( const Validator<ATTRIBUTE_T>& validator );

  void setDefaultValue ( const ATTRIBUTE_T& defaultValue );

  void readValue ( XMLReader* xmlReader );

  void readValueSpecific (
    XMLReader* xmlReader,
    double&    value );

  void readValueSpecific (
    XMLReader* xmlReader,
    int&       value );

  void readValueSpecific (
    XMLReader*   xmlReader,
    std::string& value );

  void readValueSpecific (
    XMLReader* xmlReader,
    bool&      value );

  void readValueSpecific (
    XMLReader*       xmlReader,
    utils::Vector2D& value );

  void readValueSpecific (
    XMLReader*       xmlReader,
    utils::Vector3D& value );

  void readValueSpecific (
    XMLReader*        xmlReader,
    utils::DynVector& value );

  const std::string& getName() const { return _name; };

  const ATTRIBUTE_T& getValue() const { return _value; };

  void setRead ( bool read ) { _read = read; };

  bool isRead () const { return _read; };

  /**
   * @brief Returns a documentation string about the attribute.
   */
  std::string printDocumentation() const;

private:

  static tarch::logging::Log _log;

  std::string _name;

  std::string _doc;

  bool _read;

  ATTRIBUTE_T _value;

  bool _hasDefaultValue;

  ATTRIBUTE_T _defaultValue;

  bool _hasValidation;

  Validator<ATTRIBUTE_T>* _validator;

  /**
   * @brief Sets non utils::DynVector type values.
   */
  template<typename VALUE_T>
  typename std::enable_if<
    tarch::utils::IsEqual<VALUE_T,ATTRIBUTE_T>::value &&
    not tarch::utils::IsEqual<VALUE_T,utils::DynVector>::value,
    void
  >::type set ( ATTRIBUTE_T& toSet, const VALUE_T& setter );

  /**
   * @brief Sets utils::DynVector type values by clearing and append.
   */
  template<typename VALUE_T>
  typename std::enable_if<
    tarch::utils::IsEqual<VALUE_T,ATTRIBUTE_T>::value &&
    tarch::utils::IsEqual<VALUE_T,utils::DynVector>::value,
    void
  >::type set ( ATTRIBUTE_T& toSet, const VALUE_T& setter );
};

template<typename ATTRIBUTE_T>
tarch::logging::Log XMLAttribute<ATTRIBUTE_T>:: _log ("precice::utils::XMLAttribute");

template<typename ATTRIBUTE_T>
XMLAttribute<ATTRIBUTE_T>:: XMLAttribute()
{
  assertion(false);
}


template<typename ATTRIBUTE_T>
XMLAttribute<ATTRIBUTE_T>:: XMLAttribute
(
  const std::string& name )
:
  _name(name),
  _doc(),
  _read(false),
  _value(),
  _hasDefaultValue(false),
  _defaultValue(),
  _hasValidation(false),
  _validator(NULL)
{}

template<typename ATTRIBUTE_T>
XMLAttribute<ATTRIBUTE_T>:: XMLAttribute
(
  const XMLAttribute<ATTRIBUTE_T>& rhs )
:
  _name(rhs._name),
  _doc(rhs._doc),
  _read(rhs._read),
  _value(rhs._value),
  _hasDefaultValue(rhs._hasDefaultValue),
  _defaultValue(rhs._defaultValue),
  _hasValidation(false),
  _validator(NULL)
{
  if (rhs._hasValidation) {
    assertion(rhs._validator != NULL);
    _validator = &((rhs._validator)->clone());
    _hasValidation = true;
  }
}

template<typename ATTRIBUTE_T>
void XMLAttribute<ATTRIBUTE_T>:: setDocumentation
(
  const std::string& documentation )
{
  _doc = documentation;
}

template<typename ATTRIBUTE_T>
void XMLAttribute<ATTRIBUTE_T>:: setValidator
(
  const Validator<ATTRIBUTE_T>& validator )
{
  if (_validator) {
    delete _validator;
  }
  _validator = &(validator.clone());
  _hasValidation = true;
}

//template<typename ATTRIBUTE_T>
//typename tarch::utils::EnableIf<
//  tarch::utils::IsEqual<ATTRIBUTE_T,utils::DynVector>::value,
//  void
//>::Type XMLAttribute<ATTRIBUTE_T>:: setDefaultValue
//(
//  const utils::DynVector& defaultValue  )
//{
//  preciceTrace1("setDefaultValue()", defaultValue);
//  _hasDefaultValue = true;
//  _defaultValue.clear();
//  _defaultValue.append(defaultValue);
//}

template<typename ATTRIBUTE_T>
void XMLAttribute<ATTRIBUTE_T>:: setDefaultValue
(
  const ATTRIBUTE_T& defaultValue )
{
  preciceTrace1("setDefaultValue()", defaultValue);
  _hasDefaultValue = true;
  set(_defaultValue, defaultValue);
}

template<typename ATTRIBUTE_T>
void XMLAttribute<ATTRIBUTE_T>:: readValue
(
  XMLReader* xmlReader )
{
  preciceTrace1("readValue()", _name);
  if (_read) throw "Attribute \"" + _name + "\" is defined multiple times";
  if (xmlReader->getAttributeValue(getName().c_str()) == 0) {
    if (not _hasDefaultValue){
      throw "Attribute \"" + _name + "\" missing";
    }
    set(_value, _defaultValue);
  }
  else {
    readValueSpecific(xmlReader, _value);
    if (_validator) {
      if (not _validator->validateValue(_value)){
        std::ostringstream stream;
        stream << "Invalid value \"" << _value << "\" of attribute \""
               << getName() << "\": " << _validator->getErrorMessage();
        throw stream.str();
      }
    }
  }
  preciceDebug("Read valid attribute \"" << getName() << "\" value = " << _value);
  _read = true;
}

template<typename ATTRIBUTE_T>
void XMLAttribute<ATTRIBUTE_T>:: readValueSpecific
(
  XMLReader* xmlReader,
  double&    value )
{
  value = xmlReader->getAttributeValueAsDouble (_name.c_str());
}

template<typename ATTRIBUTE_T>
void XMLAttribute<ATTRIBUTE_T>:: readValueSpecific
(
  XMLReader* xmlReader,
  int&       value )
{
  value = xmlReader->getAttributeValueAsInt(_name.c_str());
}

template<typename ATTRIBUTE_T>
void XMLAttribute<ATTRIBUTE_T>:: readValueSpecific
(
  XMLReader*   xmlReader,
  std::string& value )
{
  value = std::string(xmlReader->getAttributeValue(_name.c_str()));
}

template<typename ATTRIBUTE_T>
void XMLAttribute<ATTRIBUTE_T>:: readValueSpecific
(
  XMLReader* xmlReader,
  bool&      value )
{
  value = xmlReader->getAttributeValueAsBool(_name.c_str());
}

template<typename ATTRIBUTE_T>
void XMLAttribute<ATTRIBUTE_T>:: readValueSpecific
(
  XMLReader*       xmlReader,
  utils::Vector2D& value )
{
  value = xmlReader->getAttributeValueAsDoubleVector<2>(_name.c_str());
}

template<typename ATTRIBUTE_T>
void XMLAttribute<ATTRIBUTE_T>:: readValueSpecific
(
  XMLReader*       xmlReader,
  utils::Vector3D& value )
{
  value = xmlReader->getAttributeValueAsDoubleVector<3>(_name.c_str());
}

template<typename ATTRIBUTE_T>
void XMLAttribute<ATTRIBUTE_T>:: readValueSpecific
(
  XMLReader*        xmlReader,
  utils::DynVector& value )
{
  value.clear();
  value.append(xmlReader->getAttributeValueAsDynamicDoubleVector(_name.c_str()));
}

template<typename ATTRIBUTE_T>
std::string XMLAttribute<ATTRIBUTE_T>:: printDocumentation() const
{
  std::ostringstream doc;
  doc << _name << "=\"{" << getTypeName(_value);
  if (_hasValidation){
    doc << ":" << _validator->getDocumentation();
  }
  doc << "}";
  if (_hasDefaultValue){
    doc << "(default:'" << _defaultValue << "')";
  }
  doc << "\"";
  return doc.str();
}

template<typename ATTRIBUTE_T>
template<typename VALUE_T>
typename std::enable_if<
  tarch::utils::IsEqual<VALUE_T,ATTRIBUTE_T>::value &&
  not tarch::utils::IsEqual<VALUE_T,utils::DynVector>::value,
  void
>::type XMLAttribute<ATTRIBUTE_T>:: set
(
  ATTRIBUTE_T&   toSet,
  const VALUE_T& setter )
{
  toSet = setter;
}

template<typename ATTRIBUTE_T>
template<typename VALUE_T>
typename std::enable_if<
  tarch::utils::IsEqual<VALUE_T,ATTRIBUTE_T>::value &&
  tarch::utils::IsEqual<VALUE_T,utils::DynVector>::value,
  void
>::type XMLAttribute<ATTRIBUTE_T>:: set
(
  ATTRIBUTE_T&   toSet,
  const VALUE_T& setter )
{
  toSet.clear();
  toSet.append(setter);
}

}} // namespace precice, utils

#endif /* PRECICE_UTILS_XMLATTRIBUTE_HPP_ */
