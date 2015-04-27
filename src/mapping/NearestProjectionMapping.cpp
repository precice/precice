// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "NearestProjectionMapping.hpp"
#include "query/FindClosest.hpp"

namespace precice {
namespace mapping {

tarch::logging::Log NearestProjectionMapping::
  _log ( "precice::mapping::NearestProjectionMapping" );

NearestProjectionMapping:: NearestProjectionMapping
(
  Constraint constraint,
  int        dimensions)
:
  Mapping(constraint, dimensions),
  _weights(),
  _hasComputedMapping(false)
{
  if (constraint == CONSISTENT){
    setInputRequirement(FULL);
    setOutputRequirement(VERTEX);
  }
  else {
    assertion1(constraint == CONSERVATIVE, constraint);
    setInputRequirement(VERTEX);
    setOutputRequirement(FULL);
  }
}

void NearestProjectionMapping:: computeMapping()
{
  preciceTrace2("computeMapping()", input()->vertices().size(),
                output()->vertices().size());
  if (getConstraint() == CONSISTENT){
    preciceDebug("Compute consistent mapping");
    _weights.resize(output()->vertices().size());
    for ( size_t i=0; i < output()->vertices().size(); i++ ){
      query::FindClosest findClosest(output()->vertices()[i].getCoords());
      findClosest(*input());
      assertion(findClosest.hasFound());
      const query::ClosestElement& closest = findClosest.getClosest();
      _weights[i].clear();
      for (const query::InterpolationElement& elem : closest.interpolationElements) {
        _weights[i].push_back(elem);
      }
    }
  }
  else {
    assertion1(getConstraint() == CONSERVATIVE, getConstraint());
    preciceDebug("Compute conservative mapping");
    _weights.resize(input()->vertices().size());
    for ( size_t i=0; i < input()->vertices().size(); i++ ){
      query::FindClosest findClosest(input()->vertices()[i].getCoords());
      findClosest(*output());
      assertion(findClosest.hasFound());
      const query::ClosestElement& closest = findClosest.getClosest();
      _weights[i].clear();
      for (const query::InterpolationElement& elem : closest.interpolationElements) {
        _weights[i].push_back(elem);
      }
    }
  }
  _hasComputedMapping = true;
}

bool NearestProjectionMapping:: hasComputedMapping()
{
  return _hasComputedMapping;
}

void NearestProjectionMapping:: clear()
{
  preciceTrace("clear()");
  _weights.clear();
  _hasComputedMapping = false;
}

void NearestProjectionMapping:: map
(
  int inputDataID,
  int outputDataID )
{
  preciceTrace2("map()", inputDataID, outputDataID);
  mesh::PtrData inData = input()->data(inputDataID);
  mesh::PtrData outData = output()->data(outputDataID);
  const utils::DynVector& inValues = inData->values();
  utils::DynVector& outValues = outData->values();
  //assign(outValues) = 0.0;
  int dimensions = inData->getDimensions();
  assertion(dimensions == outData->getDimensions());

  if (getConstraint() == CONSISTENT){
    preciceDebug("Map consistent");
    assertion2(_weights.size() == output()->vertices().size(),
               _weights.size(), output()->vertices().size());
    for (size_t i=0; i < output()->vertices().size(); i++){
      InterpolationElements& elems = _weights[i];
      size_t outOffset = i * dimensions;
      for (query::InterpolationElement& elem : elems) {
        size_t inOffset = (size_t)elem.element->getID() * dimensions;
        for (int dim=0; dim < dimensions; dim++){
          assertion(outOffset + dim < (size_t)outValues.size());
          assertion(inOffset + dim < (size_t)inValues.size());
          outValues[outOffset + dim] += elem.weight * inValues[inOffset + dim];
        }
      }
    }
  }
  else {
    assertion1(getConstraint() == CONSERVATIVE, getConstraint());
    preciceDebug("Map conservative");
    assertion2(_weights.size() == input()->vertices().size(),
               _weights.size(), input()->vertices().size());
    for (size_t i=0; i < input()->vertices().size(); i++){
      size_t inOffset = i * dimensions;
      InterpolationElements& elems = _weights[i];
      for (query::InterpolationElement& elem : elems) {
        size_t outOffset = (size_t)elem.element->getID() * dimensions;
        for ( int dim=0; dim < dimensions; dim++ ){
          assertion(outOffset + dim < (size_t)outValues.size());
          assertion(inOffset + dim < (size_t)inValues.size());
          outValues[outOffset + dim] += elem.weight * inValues[inOffset + dim];
        }
      }
    }
  }
}

bool NearestProjectionMapping::doesVertexContribute(
  int vertexID)
{
  preciceTrace1("doesVertexContribute()", vertexID);
  if (getConstraint() == CONSISTENT) {
    for (size_t i=0; i < output()->vertices().size(); i++) {
      InterpolationElements& elems = _weights[i];
      for (query::InterpolationElement& elem : elems) {
        if (elem.element->getID()==vertexID && elem.weight!=0.0) {
          return true;
        }
      }
    }
  }
  else {
    assertion1(getConstraint() == CONSERVATIVE, getConstraint());
    for (size_t i=0; i < input()->vertices().size(); i++) {
      InterpolationElements& elems = _weights[i];
      for (query::InterpolationElement& elem : elems) {
        if (elem.element->getID()==vertexID){ // && elem.weight!=0.0){
          return true;
        }
      }
    }
  }
  return false;
}

}} // namespace precice, mapping
