#pragma once

#include <Eigen/Core>
#include <vector>
#include "cplscheme/SharedPointer.hpp"

namespace precice {
namespace com {
class Communication;

namespace serialize {

/// serialized representation of CouplingData
class SerializedStamples {
public:
  /**
   * @brief Serializes a given CouplingData into SerializedStamples
   *
   * @param data pointer to CouplingData to be serialized
   * @return SerializedStamples contains the serialized data
   */
  static SerializedStamples serialize(const cplscheme::PtrCouplingData data);

  /**
   * @brief Create SerializedStamples with allocated buffers according to size of CouplingData
   *
   * @param timeStamps Corresponding time stamps that will be stored in SerializedSamples
   * @param data pointer to CouplingData defining size of buffer and whether gradient data exists
   * @return SerializedStamples has allocated data buffers for serialized data
   */
  static SerializedStamples empty(Eigen::VectorXd timeStamps, const cplscheme::PtrCouplingData data);

  /**
   * @brief Deserialize data from this SerializedStamples into provided CouplingData
   *
   * @param timeStamps Corresponding time stamps for deserialized data
   * @param data pointer to CouplingData the SerializedStampes will be deserialized into
   */
  void deserializeInto(Eigen::VectorXd timeStamps, const cplscheme::PtrCouplingData data);

  /**
   * @brief const reference to serialized values. Used for sending serialized values.
   *
   * @return const Eigen::VectorXd&
   */
  const Eigen::VectorXd &values() const;

  /**
   * @brief Reference to serialized gradients. Used for storing received serialized values into.
   *
   * @return const Eigen::VectorXd&
   */
  Eigen::VectorXd &values();

  /**
   * @brief const reference to serialized gradients. Used for sending serialized gradients.
   *
   * @return const Eigen::VectorXd&
   */
  const Eigen::VectorXd &gradients() const;

  /**
   * @brief Reference to serialized gradients. Used for storing received serialized values into.
   *
   * @return const Eigen::VectorXd&
   */
  Eigen::VectorXd &gradients();

  /**
 * @brief Returns number of timeSteps
 *
 * @return int number of time steps
 */
  int nTimeSteps() const;

private:
  SerializedStamples() = default;

  // Allocates _values and gradients for size matching data and _timeSteps
  void allocate(const cplscheme::PtrCouplingData data);

  /**
   * @brief Serialize values from timeStepsStorage of data into _values
   *
   * @param data
   */
  void serializeValues(const cplscheme::PtrCouplingData data);

  void serializeValuesInitialization(const cplscheme::PtrCouplingData data);

  /**
   * @brief Serialize gradients from timeStepsStorage of data into _gradients
   *
   * @param data
   */
  void serializeGradients(const cplscheme::PtrCouplingData data);

  void serializeGradientsInitialization(const cplscheme::PtrCouplingData data);

  /**
     * @brief Deserialize _values and (if required by data) _gradients into  timeStepsStorage of data. Use provided timeStamps.
     *
     * @param timeStamps
     * @param data
     */
  void deserialize(const Eigen::VectorXd timeStamps, cplscheme::PtrCouplingData data) const;

  /// Buffer for serialized values of stamples
  Eigen::VectorXd _values;

  /// Buffer for serialized gradients of stamples
  Eigen::VectorXd _gradients;

  /// number of timesteps stored in SerializedStamples
  int _timeSteps;
};

} // namespace serialize
} // namespace com
} // namespace precice
