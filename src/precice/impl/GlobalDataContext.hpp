#pragma once

#include <string>
// #include "MappingContext.hpp"
// #include "MeshContext.hpp"
#include "mesh/SharedPointer.hpp"

namespace precice {

namespace testing {
// Forward declaration to friend the boost test struct
class DataContextFixture;
} // namespace testing
//TODO: May need a GlobalDataContextFixture

namespace impl {

/**
 * @brief Stores one Global (meshless) Data object.
 *
 * - GlobalDataContext is basically a stripped-down equivalent of DataContext. 
 *   Unlike DataContext, it has no associated mappings or meshes.
 */
class GlobalDataContext {
  friend class testing::DataContextFixture; // Make the fixture friend of this class
//TODO: May need a GlobalDataContextFixture

public:
  /**
   * @brief Get the Name of _providedData.
   *
   * @return std::string Name of _providedData.
   */
  std::string getDataName() const;

  /**
   * @brief Get the dimensions of _providedData.
   *
   * @return int Dimensions of _providedData.
   */
  int getDataDimensions() const;

  /**
   * @brief Resets provided data and (if mapping exists) fromData or toData.
   */
  void resetData();

protected:
  /**
   * @brief Construct a new DataContext without a mapping. Protected, because only ReadDataContext and WriteDataContext should use this constructor.
   *
   * @param data Data associated with this DataContext.
   * @param mesh Mesh associated with this DataContext.
   */
  GlobalDataContext(mesh::PtrGlobalData data);

  /// Unique data this context is associated with
  mesh::PtrGlobalData _providedData;

private:
  static logging::Logger _log;

};

} // namespace impl
} // namespace precice
