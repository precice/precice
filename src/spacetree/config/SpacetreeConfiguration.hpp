// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_NEWSPACETREE_SPACETREECONFIGURATION_HPP_
#define PRECICE_NEWSPACETREE_SPACETREECONFIGURATION_HPP_

#include "utils/Dimensions.hpp"
#include "utils/xml/XMLTag.hpp"
#include "spacetree/SharedPointer.hpp"
#include <string>
#include <vector>

namespace precice {
  namespace spacetree {
    class Spacetree;
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace spacetree {

/**
 * @brief Parses the xml tag of a spacetree and stores the found information.
 */
class SpacetreeConfiguration : public utils::XMLTag::Listener
{
public:

  struct ConfiguredSpacetree
  {
    std::string name;
    PtrSpacetree spacetree;
  };

  // @brief Name of the xml-tag corresponding to the SpacetreeConfiguration.
  //static const std::string& getTag();

  /**
   * @brief Constructor.
   */
  SpacetreeConfiguration ( utils::XMLTag& parent );

  /**
   * @brief Destructor, empty.
   */
  virtual ~SpacetreeConfiguration() {}

  /**
   * @brief Sets the spatial dimensionality for the spacetrees to be configured.
   */
  void setDimensions ( int dimensions );

  /**
   * Parse a TAG-tag. The parser reads the tag, validates the context, ensures
   * the corresponding compiler switch is set or not (depending on tag) and
   * returns. The argument may not be 0. If the either the validation or the
   * compiler switch check fails, all successing isValid() calls fail.
   */
  //bool parseSubtag ( utils::XMLTag::XMLReader* xmlReader );

  /**
   * @brief Returns a factory that creates the configured spacetree.
   */
//  PtrSpacetree getSpacetree () const;

  const PtrSpacetree& getSpacetree ( const std::string& name ) const;

  /**
   * @brief Returns, whether all compiler switches checked using the parseSubtag()
   *        routine were valid compared to the values specified in the
   *        configuration file.
   */
//  bool isValid() const
//  {
//    return _isValid;
//  }

  /**
   * @brief Callback function required for use of automatic configuration
   */
  virtual void xmlTagCallback ( utils::XMLTag& callingTag );

  /**
   * @brief Callback function required for use of automatic configuration.
   *
   * Is called by utils::XMLTag on automatic configuration every time an xml
   * end tag is reached.
   * @param callingTag [IN] XML tag whose end has been reached.
   * @param xmlReader  [IN] XML reader responsible for reading the tag.
   * @return True, if the callback has been successful.
   */
  virtual void xmlEndTagCallback ( utils::XMLTag& callingTag ) {}

private:

  static tarch::logging::Log _log;

  // @brief Name of the xml-tag read by this class
  const std::string TAG;
  const std::string ATTR_NAME;
  const std::string ATTR_TYPE;
  const std::string VALUE_DYNAMIC_OCTREE;
  const std::string VALUE_STATIC_OCTREE;
  const std::string VALUE_DYNAMIC_PEANOTREE2D;
  const std::string VALUE_DYNAMIC_PEANOTREE3D;

  // @brief Spatial dimensions of the spacetree to be configured.
  int _dimensions;

  std::vector<ConfiguredSpacetree> _spacetrees;

  // @brief True, if the configuration has taken place and is valid
  //bool _isValid;

  PtrSpacetree getSpacetree (
    const std::string&      type,
    const utils::DynVector& offset,
    const utils::DynVector& halflengths,
    double                  maxMeshwidth ) const;
};

}} // namespace precice, spacetree

#endif // PRECICE_NEWSPACETREE_SPACETREECONFIGURATION_HPP_
