// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_UNCOUPLEDCOUPLINGSCHEME_HPP_
#define PRECICE_CPLSCHEME_UNCOUPLEDCOUPLINGSCHEME_HPP_

#include "BaseCouplingScheme.hpp"
#include "tarch/logging/Log.h"

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace cplscheme {

class UncoupledScheme : public BaseCouplingScheme
{
public:

  /**
   * @brief Constructor.
   */
  UncoupledScheme (
    double             maxTime,
    int                maxTimesteps,
    int                validDigits,
    const std::string& participant );

  /**
   * @brief Destructor, virtual because of virtual functions.
   */
  virtual ~UncoupledScheme() {}

  /**
   * @brief Initializes the coupling scheme.
   *
   * Sets up communication to coupling partner, initializes coupling state.
   */
  virtual void initialize (
    double startTime,
    int    startTimestep );

  /**
   * @brief Empty, since this scheme has no data.
   */
  virtual void initializeData();

  /**
   * @brief Adds newly computed time. Has to be called before every advance.
   */
  virtual void addComputedTime ( double timeToAdd );

  /**
   * @brief Advances within the coupling scheme.
   */
  virtual void advance();

  /**
   * @brief Empty, since nothing to finalize.
   */
  virtual void finalize();

  /*
   * @brief returns list of all coupling partners
   */
  virtual std::vector<std::string> getCouplingPartners () const;

  virtual void sendState (
    com::Communication::SharedPointer communication,
    int                   rankReceiver );

  virtual void receiveState (
    com::Communication::SharedPointer communication,
    int                   rankSender );

  virtual std::string printCouplingState() const;

  virtual void exportState(const std::string& filenamePrefix) const {}

  virtual void importState(const std::string& filenamePrefix) {}

private:

  // @brief Logging device.
  static tarch::logging::Log _log;

  // @brief Name of participant starting the explicit coupling scheme.
  std::string _participant;
};

}} // namespace precice, cplscheme

#endif /* PRECICE_CPLSCHEME_UNCOUPLEDCOUPLINGSCHEME_HPP_ */
