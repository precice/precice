// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "Mapping.hpp"

namespace precice {
namespace mapping {

tarch::logging::Log Mapping:: _log ( "precice::mapping::Mapping" );

Mapping:: Mapping
(
  Constraint      constraint,
  int             dimensions)
:
  _constraint(constraint),
  _inputRequirement(UNDEFINED),
  _outputRequirement(UNDEFINED),
  _input(),
  _output(),
  _dimensions(dimensions)
{}

void Mapping:: setMeshes
(
  const mesh::PtrMesh& input,
  const mesh::PtrMesh& output )
{
  _input = input;
  _output = output;
}

const mesh::PtrMesh& Mapping:: getInputMesh(){
  return _input;
}

const mesh::PtrMesh& Mapping:: getOutputMesh(){
  return _output;
}

Mapping::Constraint Mapping:: getConstraint() const
{
  return _constraint;
}

Mapping::MeshRequirement Mapping:: getInputRequirement() const
{
  return _inputRequirement;
}

Mapping::MeshRequirement Mapping:: getOutputRequirement() const
{
  return _outputRequirement;
}

mesh::PtrMesh Mapping:: input()
{
  return _input;
}

mesh::PtrMesh Mapping:: output()
{
  return _output;
}

void Mapping:: setInputRequirement
(
  MeshRequirement requirement )
{
  _inputRequirement = requirement;
}

void Mapping:: setOutputRequirement
(
  MeshRequirement requirement )
{
  _outputRequirement = requirement;
}

bool Mapping:: doesVertexContribute(
  int vertexID)
{
  return true;
}

int Mapping:: getDimensions(){
  return _dimensions;
}

}} // namespace precice, mapping

