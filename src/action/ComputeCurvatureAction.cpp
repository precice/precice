// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#include "ComputeCurvatureAction.hpp"
#include "utils/Globals.hpp"
#include "utils/Dimensions.hpp"
#include "mesh/Mesh.hpp"
#include "mesh/Triangle.hpp"
#include "mesh/Edge.hpp"
#include "mesh/Vertex.hpp"
#include "utils/Globals.hpp"

namespace precice {
namespace action {

tarch::logging::Log ComputeCurvatureAction::
  _log ( "precice::action::ComputeCurvatureAction" );

ComputeCurvatureAction:: ComputeCurvatureAction
(
  Timing               timing,
  int                  dataID,
  const mesh::PtrMesh& mesh )
:
  Action(timing, mesh),
  _data(mesh->data(dataID))
{}

void ComputeCurvatureAction:: performAction
(
  double time,
  double dt,
  double computedPartFullDt,
  double fullDt )
{
  preciceTrace ( "performAction()" );
  utils::DynVector& dataValues = _data->values();

  if ( getMesh()->getDimensions() == 2 ){
    assign(dataValues) = 0.0;
    utils::Vector2D tangent;
    foreach ( mesh::Edge & edge, getMesh()->edges() ){
      mesh::Vertex& v0 = edge.vertex(0);
      mesh::Vertex& v1 = edge.vertex(1);
      tangent = v1.getCoords();
      tangent -= v0.getCoords();
      tangent /= norm2(tangent);
      for(int d=0; d < 2; d++) {
        dataValues[v0.getID()*2 + d] += tangent[d];
        dataValues[v1.getID()*2 + d] -= tangent[d];
      }
    }
  }
  else {
    assertion1 ( getMesh()->getDimensions() == 3, getMesh()->getDimensions() );
    for (int i=0; i < dataValues.size(); i++) {
      dataValues[i] = 0.0;
    }
    utils::Vector3D normal;
    utils::Vector3D edge;
    utils::Vector3D contribution;

    foreach ( mesh::Triangle& tri, getMesh()->triangles() ){
      normal = tri.getNormal();
      for (int i=0; i < 3; i++) {
        mesh::Vertex& v0 = tri.vertex(i);
        mesh::Vertex& v1 = tri.vertex((i+1) % 3);
        edge = v1.getCoords();
        edge -= v0.getCoords();
        cross (edge, normal, contribution);
        for ( int d=0; d < 3; d++ ) {
          dataValues[v0.getID() * 3 + d] -= 0.25 * contribution[d];
          dataValues[v1.getID() * 3 + d] -= 0.25 * contribution[d];
        }
      }
    }
  }
}

}} // namespace precice, action
