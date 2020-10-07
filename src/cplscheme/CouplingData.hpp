#pragma once

#include <Eigen/Core>
#include "mesh/Data.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/assertion.hpp"
#include "utils/EigenHelperFunctions.hpp"

namespace precice {
namespace cplscheme {

struct CouplingData {  // @todo: should be a class from a design standpoint. See https://github.com/precice/precice/pull/865#discussion_r495825098
  using DataMatrix = Eigen::MatrixXd;

  /// Returns a reference to the data values.
  Eigen::VectorXd &values()
  {
    PRECICE_ASSERT(data != nullptr);
    return data->values();
  }

  /// Returns a const reference to the data values.
  const Eigen::VectorXd &values() const
  {
    PRECICE_ASSERT(data != nullptr);
    return data->values();
  }

  void updateOldValues()
  {
    oldValues.col(0) = this->values();
    // For extrapolation, treat the initial value as old time windows value
    utils::shiftSetFirst(this->oldValues, this->values());
  }

  void extrapolateData(int order, int timeWindows)
  {
    if ((order == 1) || (timeWindows == 2)) { //timesteps is increased before extrapolate is called
      if((order != 1) && (timeWindows == 2))
        PRECICE_ASSERT(order == 2);
      // PRECICE_INFO("Performing first order extrapolation");
      PRECICE_ASSERT(this->oldValues.cols() > 1);
      this->values() *= 2.0;              // = 2*x^t
      this->values() -= this->oldValues.col(1); // = 2*x^t - x^(t-1)
      utils::shiftSetFirst(this->oldValues, this->values());
    } else if (order == 2) {
      // PRECICE_INFO("Performing second order extrapolation");
      PRECICE_ASSERT(this->oldValues.cols() > 2);
      auto valuesOld1 = this->oldValues.col(1);
      auto valuesOld2 = this->oldValues.col(2);

      this->values() *= 2.5;              // = 2.5 x^t
      this->values() -= valuesOld1 * 2.0; // = 2.5x^t - 2x^(t-1)
      this->values() += valuesOld2 * 0.5; // = 2.5x^t - 2x^(t-1) + 0.5x^(t-2)
      utils::shiftSetFirst(this->oldValues, this->values());
    } else {
      PRECICE_ASSERT(false, "Extrapolation order is invalid.");
    }
  }

  /// Data values of previous iteration (1st col) and previous time windows.
  DataMatrix oldValues;

  mesh::PtrData data;

  mesh::PtrMesh mesh;

  ///  True, if the data values if this CouplingData requires to be initialized by a participant.
  bool requiresInitialization;

  int getDimensions()
  {
    PRECICE_ASSERT(data != nullptr);
    return data->getDimensions();
  }

  /**
   * @brief Default constructor, not to be used!
   *
   * Necessary when compiler creates template code for std::map::operator[].
   */
  CouplingData()
  {
    PRECICE_ASSERT(false);
  }

  CouplingData(
      mesh::PtrData data,
      mesh::PtrMesh mesh,
      bool          requiresInitialization)
      : data(data),
        mesh(mesh),
        requiresInitialization(requiresInitialization)
  {
    PRECICE_ASSERT(data != nullptr);
    PRECICE_ASSERT(mesh != nullptr);
    PRECICE_ASSERT(mesh.use_count() > 0);
  }
};

} // namespace cplscheme
} // namespace precice
