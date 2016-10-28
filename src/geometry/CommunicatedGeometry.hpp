#pragma once

#include "Geometry.hpp"
#include "m2n/M2N.hpp"
#include "mapping/SharedPointer.hpp"
#include "utils/MasterSlave.hpp"
#include "logging/Logger.hpp"
#include "impl/Decomposition.hpp"
#include "impl/SharedPointer.hpp"
#include <string>
#include <map>

namespace precice {
namespace geometry {

/**
 * @brief Creates the mesh by copying another remote mesh via communication. In case of parallel
 * solvers using a master-slave approach the meshes are gathered or scattered. In case of scattering,
 * globals meshes are sent first to each slave and filtered then depending on both mappings defined.
 */
class CommunicatedGeometry : public Geometry
{
public:

  CommunicatedGeometry (
    const Eigen::VectorXd&  offset,
    const std::string&      accessor,
    const std::string&      provider,
    impl::PtrDecomposition  decomposition);

  virtual ~CommunicatedGeometry() {}

  void addReceiver (
    const std::string& receiver,
    m2n::M2N::SharedPointer m2n );

  /*
   * @bried Prepare geometry for creation, i.e. communicate mesh
   */
  void prepare ( mesh::Mesh& seed );

protected:

  /**
   * Is called from Geometry and sends the mesh if the accessor is provider,
   * receives if the accessor is contained in receivers, fails otherwise.
   */
  void specializedCreate ( mesh::Mesh& seed );

private:

  void sendMesh(
    mesh::Mesh& seed);

  void receiveMesh(
    mesh::Mesh& seed);

  /// Logging device.
  static logging::Logger _log;

  std::string _accessorName;

  std::string _providerName;

  std::map<std::string,m2n::M2N::SharedPointer> _receivers;

  impl::PtrDecomposition _decomposition;
};

}} // namespace precice, geometry
