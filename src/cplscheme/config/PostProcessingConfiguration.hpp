// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_POSTPROCESSINGCONFIGURATION_HPP_
#define PRECICE_CPLSCHEME_POSTPROCESSINGCONFIGURATION_HPP_

#include "cplscheme/impl/SharedPointer.hpp"
#include "cplscheme/impl/PostProcessing.hpp"
#include "cplscheme/impl/MVQNPostProcessing.hpp"
#include "precice/config/SharedPointer.hpp"
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
    * @brief Returns a pointer to the PostProcessingConfig object for the coarse model optimization method
    */
   PtrPostProcessingConfiguration getCoarseModelOptimizationConfig();

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

   void setIsAddManifoldMappingTagAllowed(bool b)
   {
     _isAddManifoldMappingTagAllowed = b;
   }

private:

   static tarch::logging::Log _log;

   const std::string TAG;
   const std::string TAG_RELAX;
   const std::string TAG_INIT_RELAX;
   const std::string TAG_MAX_USED_ITERATIONS;
   const std::string TAG_TIMESTEPS_REUSED;
   const std::string TAG_DATA;
   const std::string TAG_FILTER;
   const std::string TAG_ESTIMATEJACOBIAN;
   const std::string TAG_PRECONDITIONER;
   const std::string TAG_IMVJRESTART;

   const std::string ATTR_NAME;
   const std::string ATTR_MESH;
   const std::string ATTR_SCALING;
   const std::string ATTR_VALUE;
   const std::string ATTR_ENFORCE;
   const std::string ATTR_SINGULARITYLIMIT;
   const std::string ATTR_TYPE;
   const std::string ATTR_BUILDJACOBIAN;
   const std::string ATTR_IMVJCHUNKSIZE;
   const std::string ATTR_RSLS_REUSEDTSTEPS;
   const std::string ATTR_RSSVD_TRUNCATIONEPS;
   const std::string ATTR_PRECOND_NONCONST_TIMESTEPS;

   const std::string VALUE_CONSTANT;
   const std::string VALUE_AITKEN;
   const std::string VALUE_HIERARCHICAL_AITKEN;
   const std::string VALUE_IQNILS;
   const std::string VALUE_MVQN;
   const std::string VALUE_ManifoldMapping;
   const std::string VALUE_BROYDEN;
   const std::string VALUE_QR1FILTER;
   const std::string VALUE_QR1_ABSFILTER;
   const std::string VALUE_QR2FILTER;
   const std::string VALUE_CONSTANT_PRECONDITIONER;
   const std::string VALUE_VALUE_PRECONDITIONER;
   const std::string VALUE_RESIDUAL_PRECONDITIONER;
   const std::string VALUE_RESIDUAL_SUM_PRECONDITIONER;
   const std::string VALUE_LS_RESTART;
   const std::string VALUE_ZERO_RESTART;
   const std::string VALUE_SVD_RESTART;
   const std::string VALUE_NO_RESTART;

   //bool _isValid;

   const mesh::PtrMeshConfiguration _meshConfig;

   std::string _meshName;

   // post processing method
   impl::PtrPostProcessing _postProcessing;

   // recursive definition of post processings for multi level methods (i.e., manifold mapping)
   PtrPostProcessingConfiguration _coarseModelOptimizationConfig;

   std::vector<std::string> _neededMeshes;

   impl::PtrPreconditioner _preconditioner;

   struct ConfigurationData
   {
      std::vector<int> dataIDs;
      std::map<int,double> scalings;
      std::string type;
      double relaxationFactor;
      bool forceInitialRelaxation;
      int maxIterationsUsed;
      int timestepsReused;
      int filter;
      int imvjRestartType;
      int imvjChunkSize;
      int imvjRSLS_reustedTimesteps;
      int precond_nbNonConstTSteps;
      double singularityLimit;
      double imvjRSSVD_truncationEps;
      bool estimateJacobian;
      bool alwaysBuildJacobian;
      std::string preconditionerType;

      ConfigurationData ()
      :
         dataIDs (),
         scalings(),
         type ( "" ),
         relaxationFactor ( 0.0 ),
         forceInitialRelaxation( false ),
         maxIterationsUsed ( 0 ),
         timestepsReused ( 0 ),
         filter ( impl::PostProcessing::NOFILTER ),
         imvjRestartType( 0 ), // NO-RESTART
         imvjChunkSize ( 0 ),
         imvjRSLS_reustedTimesteps( 0 ),
         precond_nbNonConstTSteps( -1),
         singularityLimit ( 0.0 ),
         imvjRSSVD_truncationEps( 0.0 ),
         estimateJacobian ( false ),
         alwaysBuildJacobian( false ),
         preconditionerType("")
      {}

   } _config;

   bool _isAddManifoldMappingTagAllowed;


   void addTypeSpecificSubtags ( utils::XMLTag& tag );
};

}} // namespace precice, cplscheme

#endif /* PRECICE_CPLSCHEME_POSTPROCESSINGCONFIGURATION_HPP_ */
