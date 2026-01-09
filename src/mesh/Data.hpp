#pragma once

#include <Eigen/Core>
#include <initializer_list>
#include <stddef.h>
#include <string>

#include "SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "precice/impl/Types.hpp"
#include "time/Sample.hpp"
#include "time/Time.hpp"
#include "time/Waveform.hpp"

namespace precice::mesh {
class Mesh;
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice::mesh {

/**
 * @brief Describes a set of data values belonging to the vertices of a mesh.
 */
class Data {
public:
  // @brief Data dimensions type (scalar/vector)
  enum typeName {
    SCALAR,
    VECTOR
  };

  // @brief Possible types of data values.
  //  enum DataType {
  //    TYPE_UNDEFINED,
  //    TYPE_DOUBLE,
  //    TYPE_VECTOR
  //  };

  // @brief Name of an undefined data type.
  // static const std::string TYPE_NAME_UNDEFINED;
  // @brief Name of a double data type.
  // static const std::string TYPE_NAME_DOUBLE;
  // @brief Name of a vector data type.
  // static const std::string TYPE_NAME_VECTOR;

  /**
   * @brief Constructor
   */
  Data(
      std::string                        name,
      DataID                             id,
      int                                dimension,
      int                                spatialDimensions = -1,
      int                                waveformDegree    = time::Time::DEFAULT_WAVEFORM_DEGREE,
      std::vector<std::optional<double>> lowerBound        = std::vector<std::optional<double>>(3),
      std::vector<std::optional<double>> upperBound        = std::vector<std::optional<double>>(3));

  /// Returns a reference to the data values.
  Eigen::VectorXd &values();

  /// Returns a const reference to the data values.
  const Eigen::VectorXd &values() const;

  /// Returns a const reference to the gradient data values.
  const Eigen::MatrixXd &gradients() const;

  /// Returns a const reference to the _sample.
  const time::Sample &sample() const;

  /**
   * @brief Samples _waveform at given time
   *
   * @param time Time where the sampling happens.
   * @return Value of _waveform at time \ref time.
   */
  time::SampleResult sampleAtTime(double time) const;

  /**
   * @brief get degree of _waveform.
   *
   * @return int degree of _waveform
   */
  int getWaveformDegree() const;

  /// Returns a reference to the waveform
  time::Waveform &waveform();

  std::vector<std::optional<double>> getLowerBound() const;

  std::vector<std::optional<double>> getUpperBound() const;

  void moveToNextWindow();

  /// Returns a the stamples from the waveform
  auto stamples() const
  {
    return _waveform.stamples();
  }

  /// Add sample at given time to the waveform
  void setSampleAtTime(double time, const time::Sample &sample);

  /// Set _sample
  void setGlobalSample(const time::Sample &sample); // @todo try to remove this function

  /// Creates an empty sample at given time
  void emplaceSampleAtTime(double time);

  /// Creates a sample at given time with given values
  void emplaceSampleAtTime(double time, std::initializer_list<double> values);

  /// Creates a sample at given time with given values and gradients
  void emplaceSampleAtTime(double time, std::initializer_list<double> values, std::initializer_list<double> gradients);

  /// Returns the name of the data set, as set in the config file.
  const std::string &getName() const;

  /// Returns the ID of the data set (supposed to be unique).
  DataID getID() const;

  /// Returns if the data contains gradient data
  bool hasGradient() const;

  /// Returns if there are sample of this data
  bool hasSamples() const;

  /// Set the additional requirement of gradient data
  void requireDataGradient();

  /// Returns the mesh dimension (i.e., number of rows) of one gradient data value .
  int getSpatialDimensions() const;

  /// Returns the dimension (i.e., number of components) of one data value (i.e number of columns of one gradient data value).
  int getDimensions() const;

  /**
   * @brief Allocates memory for the data values and corresponding gradient values.
   *
   * @param expectedCount expected number of values count (i.e. number of mesh vertices)
   */
  void allocateValues(int expectedCount);

private:
  logging::Logger _log{"mesh::Data"};

  /// Sample storage of this Data.
  time::Waveform _waveform;

  /// Lower bound for data values.
  std::vector<std::optional<double>> _lowerBound;

  /// Upper bound for data values.
  std::vector<std::optional<double>> _upperBound;

  /// Name of the data set.
  std::string _name;

  /// ID of the data set (supposed to be unique).
  DataID _id;

  /// Dimensionality of one data value.
  int _dimensions;

  /// Spatial Dimension of one element -> number of rows (only 2, 3 allowed for 2D, 3D).
  int _spatialDimensions;

  /// Whether gradient data is available or not
  bool _hasGradient = false;

  time::Sample _sample;
};

} // namespace precice::mesh
