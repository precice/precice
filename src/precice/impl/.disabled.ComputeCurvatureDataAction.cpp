// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ComputeCurvatureDataAction.hpp"

#include "utils/Globals.hpp"
#include "utils/Dimensions.hpp"

#include "mesh/Edge.hpp"
#include "mesh/Vertex.hpp"

namespace precice {
namespace impl {

tarch::logging::Log ComputeCurvatureDataAction::
  _log ( "precice::impl::ComputeCurvatureDataAction" );

ComputeCurvatureDataAction:: ComputeCurvatureDataAction
(
  TimingConstants       timing,
  int                   dataID,
  const mesh::PtrMesh & mesh )
:
  AbstractDataAction ( timing ),
  _data ( mesh->data(dataID) ),
  _mesh ( mesh )
{}

void ComputeCurvatureDataAction:: performAction
(
  double dt,
  double computedPartFullDt,
  double fullDt )
{
  preciceTrace ( "performAction()" );

  preciceDebug ( "Computing curvature" );


#ifdef Dim2
  mesh::Data::Values& dataValues = _data->values();

  for (int i=0; i<dataValues.size(); i++)
    dataValues(i) = 0.0;

  foreach( mesh::Edge & edge, _mesh->edges() ) {
    mesh::Vertex& v0 = edge.vertex(0);
    mesh::Vertex& v1 = edge.vertex(1);
    ::precice::utils::Vector tangent = v1.getCoords() - v0.getCoords();
    tangent /= norm2(tangent);

    for(int i=0; i<PRECICE_DIMENSIONS; i++) {
      dataValues(v0.getID()*PRECICE_DIMENSIONS+i) += tangent(i);
      dataValues(v1.getID()*PRECICE_DIMENSIONS+i) -= tangent(i);
    }
  }
#elif Dim3
  preciceError("ComputeCurvatureDataAction::performAction()",
                     "Curvature computation for 3D not implemented yet");
#endif
}

}} // namespace precice, impl
