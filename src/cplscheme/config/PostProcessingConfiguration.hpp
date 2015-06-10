// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_POSTPROCESSINGCONFIGURATION_HPP_
#define PRECICE_CPLSCHEME_POSTPROCESSINGCONFIGURATION_HPP_

#include "cplscheme/impl/SharedPointer.hpp"
#include "mesh/SharedPointer.hpp"
#include "utils/xml/XMLTag.hpp"
#include "tarch/logging/Log.h"

namespace precice {
namespace cplscheme {

class PostProcessingConfiguration : public utils::XMLTag::Listener
{
public:

   // @brief Name of the XML tag holding post-processing information.
   //static const std::string& getTag();

   /**
    * @brief Constructor.
    */
   PostProcessingConfiguration (
     const mesh::PtrMeshConfiguration& meshConfig);

   /**
    * @brief Parses the XML information in xmlReader to a configuration.
    */
   //bool parseSubtag ( tarch::irr::io::IrrXMLReader * xmlReader );

   /**
    * @brief Returns true, if configuration has taken place and is valid.
    */
   //bool isValid() const;

   /**
    * @brief Returns the configured (to be checked by isValid()) coupling scheme.
    */
   impl::PtrPostProcessing getPostProcessing();

   /**
    * @brief Callback method required when using utils::XMLTag.
    */
   virtual void xmlTagCallback ( utils::XMLTag& callingTag );

   /**
    * @brief Callback method required when using utils::XMLTag.
    */
   virtual void xmlEndTagCallback ( utils::XMLTag& callingTag );

   /**
    * @brief Removes configured post-processing.
    */
   void clear();

   /**
   * @brief Connect tags.
   */
   void connectTags( utils::XMLTag& tag );

   std::vector<std::string>& getNeededMeshes(){
     return _neededMeshes;
   }

private:

   static tarch::logging::Log _log;

   const std::string TAG;
   const std::string TAG_RELAX;
   const std::string TAG_INIT_RELAX;
   const std::string TAG_MAX_USED_ITERATIONS;
   const std::string TAG_TIMESTEPS_REUSED;
   const std::string TAG_SINGULARITY_LIMIT;
   const std::string TAG_DATA;

   const std::string ATTR_NAME;
   const std::string ATTR_MESH;
   const std::string ATTR_SCALING;
   const std::string ATTR_VALUE;

   const std::string VALUE_CONSTANT;
   const std::string VALUE_AITKEN;
   const std::string VALUE_HIERARCHICAL_AITKEN;
   const std::string VALUE_IQNILS;
   const std::string VALUE_MVQN;
   const std::string VALUE_BROYDEN;

   //bool _isValid;

   const mesh::PtrMeshConfiguration _meshConfig;

   std::string _meshName;

   impl::PtrPostProcessing _postProcessing;

   std::vector<std::string> _neededMeshes;

   struct ConfigurationData
   {
      std::vector<int> dataIDs;
      std::map<int,double> scalings;
      std::string type;
      double relaxationFactor;
      int maxIterationsUsed;
      int timestepsReused;
      double singularityLimit;

      ConfigurationData ()
      :
         dataIDs (),
         scalings(),
         type ( "" ),
         relaxationFactor ( 0.0 ),
         maxIterationsUsed ( 0 ),
         timestepsReused ( 0 ),
         singularityLimit ( 0.0 )
      {}

   } _config;


   void addTypeSpecificSubtags ( utils::XMLTag& tag );
};

}} // namespace precice, cplscheme

#endif /* PRECICE_CPLSCHEME_POSTPROCESSINGCONFIGURATION_HPP_ */
