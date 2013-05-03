// Copyright (C) 2011 Technische Universitaet Muenchen
// This file is part of the preCICE project. For conditions of distribution and
// use, please see the license notice at http://www5.in.tum.de/wiki/index.php/PreCICE_License
#ifndef PRECICE_CPLSCHEME_COUPLINGSCHEMECONFIGURATION_HPP_
#define PRECICE_CPLSCHEME_COUPLINGSCHEMECONFIGURATION_HPP_

#include "cplscheme/CouplingScheme.hpp"
#include "cplscheme/ImplicitCouplingScheme.hpp"
#include "cplscheme/SharedPointer.hpp"
#include "cplscheme/Constants.hpp"
#include "mesh/SharedPointer.hpp"
#include "com/SharedPointer.hpp"
#include "precice/config/SharedPointer.hpp"
#include "utils/xml/XMLTag.hpp"
#include "tarch/logging/Log.h"
#include <vector>
#include <string>
#include "boost/tuple/tuple.hpp"
#include "precice/impl/MeshContext.hpp"

namespace precice {
  namespace cplscheme {
    namespace tests {
      class ImplicitCouplingSchemeTest;
    }
  }
}

// ----------------------------------------------------------- CLASS DEFINITION

namespace precice {
namespace cplscheme {

/**
 * @brief Configuration for a coupling scheme, parser for xml coupling tag.
 */
class CouplingSchemeConfiguration : public utils::XMLTag::Listener
{
public:

   /**
    * @brief Constructor.
    */
   CouplingSchemeConfiguration (
     utils::XMLTag&                            parent,
     const mesh::PtrMeshConfiguration&         meshConfig,
     const com::PtrCommunicationConfiguration& comConfig );

   /**
    * @brief Parses the XML information in xmlReader to a configuration.
    */
   //bool parseSubtag ( utils::XMLTag::XMLReader* xmlReader );

   /**
    * @brief Returns true, if configuration has taken place and is valid.
    */
   //bool isValid () const;

   bool hasCouplingScheme ( const std::string& participantName ) const;

   /**
    * @brief Returns the configured (to be checked by isValid()) coupling scheme.
    */
   const PtrCouplingScheme& getCouplingScheme ( const std::string& participantName ) const;

   /**
    * @brief Returns the name of one dataset exchanged in the coupling scheme.
    */
   const std::string& getDataToExchange ( int index ) const;

   /**
    * @brief Callback method required when using utils::XMLTag.
    */
   virtual void xmlTagCallback ( utils::XMLTag& callingTag );

   /**
    * @brief Callback method required when using utils::XMLTag.
    */
   virtual void xmlEndTagCallback ( utils::XMLTag& callingTag );

   /**
    * @brief Adds a manually configured coupling scheme for a participant.
    */
   void addCouplingScheme (
     PtrCouplingScheme  cplScheme,
     const std::string& participantName );

private:

   // @brief Logging device.
   static tarch::logging::Log _log;

   const std::string TAG;
   const std::string TAG_PARTICIPANTS;
   const std::string TAG_EXCHANGE;
   const std::string TAG_MAX_TIME;
   const std::string TAG_MAX_TIMESTEPS;
   const std::string TAG_TIMESTEP_LENGTH;
   const std::string TAG_ABS_CONV_MEASURE;
   const std::string TAG_REL_CONV_MEASURE;
   const std::string TAG_RES_REL_CONV_MEASURE;
   const std::string TAG_MIN_ITER_CONV_MEASURE;
   const std::string TAG_MAX_ITERATIONS;
   const std::string TAG_CHECKPOINT;
   const std::string TAG_EXTRAPOLATION;

   const std::string ATTR_DATA;
   const std::string ATTR_MESH;
   const std::string ATTR_PARTICIPANT;
   const std::string ATTR_INITIALIZE;
   const std::string ATTR_TYPE;
   const std::string ATTR_FIRST;
   const std::string ATTR_SECOND;
   const std::string ATTR_VALUE;
   const std::string ATTR_VALID_DIGITS;
   const std::string ATTR_METHOD;
   const std::string ATTR_LIMIT;
   const std::string ATTR_MIN_ITERATIONS;
   const std::string ATTR_NAME;
   const std::string ATTR_TIMESTEP_INTERVAL;
   const std::string ATTR_FROM;
   const std::string ATTR_SUFFICES;

   const std::string VALUE_EXPLICIT;
   const std::string VALUE_IMPLICIT;
   const std::string VALUE_UNCOUPLED;
   const std::string VALUE_FIXED;
   const std::string VALUE_FIRST_PARTICIPANT;

   struct Config
   {
      std::string type;
      std::string name;
      int checkpointTimestepInterval;
      std::string participant;
      std::string secondParticipant;
      double maxTime;
      int maxTimesteps;
      double timestepLength;
      int validDigits;
      constants::TimesteppingMethod dtMethod;
      // @brief Tuples of exchange data, mesh, and participant name.
      typedef boost::tuple<mesh::PtrData,std::string,bool> Exchange;
      std::vector<Exchange> exchanges;
      // @brief Tuples of data ID, mesh ID, and convergence measure.
      std::vector<boost::tuple<int,bool,impl::PtrConvergenceMeasure> > convMeasures;
      int maxIterations;
      //impl::PtrPostProcessing postProcessing;
      int extrapolationOrder;

      Config()
      :
         type ( "" ),
         name ( "" ),
         checkpointTimestepInterval ( -1 ),
         participant ( "" ),
         secondParticipant ( "" ),
         maxTime ( CouplingScheme::UNDEFINED_TIME ),
         maxTimesteps ( CouplingScheme::UNDEFINED_TIMESTEPS ),
         timestepLength ( CouplingScheme::UNDEFINED_TIMESTEP_LENGTH ),
         validDigits ( 16 ),
         dtMethod ( constants::FIXED_DT ),
         exchanges (),
         convMeasures (),
         maxIterations ( -1 ),
         extrapolationOrder ( 0 )
      {}

   } _config;


   mesh::PtrMeshConfiguration _meshConfig;

   com::PtrCommunicationConfiguration _comConfig;

   PtrPostProcessingConfiguration _postProcConfig;

   bool _isValid;

   std::map<std::string,PtrCouplingScheme> _couplingSchemes;

   void addTypespecifcSubtags (
      const std::string& type,
      utils::XMLTag&     tag );

   void addTagCheckpoint ( utils::XMLTag& tag );

   void addTransientLimitTags ( utils::XMLTag& tag );

   void addTagParticipants ( utils::XMLTag& tag );

   void addTagExchange ( utils::XMLTag& tag );

   void addTagAbsoluteConvergenceMeasure ( utils::XMLTag& tag );

   void addTagRelativeConvergenceMeasure ( utils::XMLTag& tag );

   void addTagResidualRelativeConvergenceMeasure ( utils::XMLTag& tag );

   void addTagMinIterationConvergenceMeasure ( utils::XMLTag& tag );

   void addBaseAttributesTagConvergenceMeasure ( utils::XMLTag& tag );

   void addTagMaxIterations ( utils::XMLTag& tag );

   void addTagExtrapolation ( utils::XMLTag& tag );

   void addTagPostProcessing ( utils::XMLTag& tag );

   void addAbsoluteConvergenceMeasure (
      const std::string & dataName,
      const std::string & meshName,
      double              limit,
      bool                suffices );

   void addRelativeConvergenceMeasure (
      const std::string & dataName,
      const std::string & meshName,
      double              limit,
      bool                suffices );

   void addResidualRelativeConvergenceMeasure (
      const std::string & dataName,
      const std::string & meshName,
      double              limit,
      bool                suffices );

   void addMinIterationConvergenceMeasure (
      const std::string & dataName,
      const std::string & meshName,
      int                 minIterations,
      bool                suffices );

   mesh::PtrData getData (
      const std::string & dataName,
      const std::string & meshName ) const;

   PtrCouplingScheme createExplicitCouplingScheme (
     const std::string & accessor ) const;

   PtrCouplingScheme createImplicitCouplingScheme (
     const std::string & accessor ) const;

   /*
    * @brief returns name of the actual scheme holder (i.e. server name)
   */
   std::string determineCouplingSchemeHolder (
     const std::string & accessorName ) const;

   constants::TimesteppingMethod getTimesteppingMethod (
     const std::string& method ) const;

   /**
    * @brief Adds configured exchange data to be sent or received to scheme.
    */
   void addDataToBeExchanged(
     CouplingScheme&    scheme,
     const std::string& accessor) const;

   friend class tests::ImplicitCouplingSchemeTest; // For whitebox tests
};

}} // namespace precice, cplscheme

#endif /* COUPLINGSCHEMECONFIGURATION_HPP_ */
