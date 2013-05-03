// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CONTAINER_SHAREDPOINTER_HPP_
#define PRECICE_CONTAINER_SHAREDPOINTER_HPP_

#include "boost/smart_ptr.hpp"

namespace precice {
namespace mesh {

class Data;
typedef boost::shared_ptr<Data>              PtrData;

class Group;
typedef boost::shared_ptr<Group>             PtrGroup;

class Mesh;
typedef boost::shared_ptr<Mesh>              PtrMesh;

class DataConfiguration;
typedef boost::shared_ptr<DataConfiguration> PtrDataConfiguration;

class MeshConfiguration;
typedef boost::shared_ptr<MeshConfiguration> PtrMeshConfiguration;
}} // namespace precice, mesh

#endif /* PRECICE_CONTAINER_SHAREDPOINTER_HPP_ */
