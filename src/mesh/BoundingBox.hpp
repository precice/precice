#pragma once
#include "logging/Logger.hpp"
#include "mesh/Vertex.hpp"
#include <vector>

namespace precice {
namespace mesh {

/// Bounding box class for mesh and partions
class BoundingBox {

public:
  /***
  * 
  * @brief Constructor.
  * 
  * @param[in] bounds Vector of minimum and maximum points in each dimension
  * @param[in] safetyFactor Factor which enlarges the bounding box
  * 
  */
  BoundingBox(std::vector<double> bounds);
  
  /***
  * 
  * @brief Constructor.
  * 
  * @param[in] dimension Dimension of the bounding box
  * 
  */
  BoundingBox(int dimension);

  /// Default constructor
  BoundingBox();

  /// Copy Constructor
  BoundingBox(const BoundingBox &bb);

  /// Destructor
  ~BoundingBox();

  /// Comparison Operator
  bool operator==(const BoundingBox& otherBB) const;

  /// Bounding box factory function
  static BoundingBox createFromData(std::vector<double> bounds);

  /***
  * 
  * @brief Set bounds manually
  * 
  * @param[in] dimension Dimension in which the bounds are given
  * @param[in] min Minimum bound
  * @param[in] max Maximum bound
  * 
  */
  void setBounds(int dimension, double min, double max);

  /// Setter of minimum bound in given direction
  void setMin(int dimension, double min);

  /// Setter of maximum bound in given direction
  void setMax(int dimension, double max);

  /// Setter safety factor
  void setSafetyFactor(double safetyFactor);

  /// Merges the bounding box with given bounding box, also enlage by _safetyFactor
  bool mergeBoundingBoxes(const BoundingBox &otherBB);

  /// Expand bounding box
  void expandTo(const Vertex& vertices);

  /// Checks if vertex in contained in _bb
  bool isVertexInBB(const Vertex &vertex);

  /// Checks whether two bounding boxes are overlapping
  bool overlapping(const BoundingBox &otherBB);

  /// Return the value for given dimension and bound type
  double getData(int dimension, int type) const;

  /// Return data as pointer, may require to use getSize()
  const double* data() const;

  /// Return data as std::vector
  const std::vector<double> dataVector() const;

  /// Getter dimension of the bounding box
  int getDimension() const;

  /// Getter of size of the bound container
  int getSize() const;

  /// Whether the bounds container is empty or not
  bool empty();

  /// Output operator for easy logging
  friend std::ostream &operator<<(std::ostream &out, const BoundingBox &bb);

private:
  logging::Logger _log{"mesh::BoundingBox"};

  /// Number of dimensions
  int    _dimensions;

  /// Safety factor to enlarge the bounding box
  double _safetyFactor{1.0};

  /// Whether this bounding box is prepared or not
  bool _prepared{false};

  /// Container of min and max points in each dimension
  std::vector<double> _bounds;
};

} // namespace mesh
} // namespace precice
