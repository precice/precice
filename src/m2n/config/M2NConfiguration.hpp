#pragma once

#include "m2n/SharedPointer.hpp"
#include "logging/Logger.hpp"
#include "xml/XMLTag.hpp"

#include <tuple>

#include <memory>
#include <string>
#include <vector>

namespace precice {
namespace m2n {

/**
 * @brief Configuration for communication channels between solvers.
 */
class M2NConfiguration : public xml::XMLTag::Listener
{
public:
   using SharedPointer = std::shared_ptr<M2NConfiguration>;

   typedef std::tuple<m2n::PtrM2N,std::string,std::string> M2NTuple;

public:
   
   M2NConfiguration ( xml::XMLTag& parent );

   virtual ~M2NConfiguration() {}

   /**
    * @brief Returns the communication object for the given user names.
    *
    * Exits with an error message, when no object is configured for the given
    * user names.
    */
   m2n::PtrM2N getM2N (
      const std::string& from,
      const std::string& to );

   /**
    * @brief Returns all configured communication objects.
    */
   std::vector<M2NTuple>& m2ns()
   {
      return _m2ns;
   }

   virtual void xmlTagCallback ( xml::XMLTag& callingTag );

   virtual void xmlEndTagCallback ( xml::XMLTag& callingTag ) {}

private:

   static logging::Logger _log;

   const std::string TAG;
   const std::string ATTR_TYPE;
   const std::string ATTR_DISTRIBUTION_TYPE;
   const std::string ATTR_FROM;
   const std::string ATTR_TO;
   const std::string ATTR_PORT;
   const std::string ATTR_NETWORK;
   const std::string ATTR_EXCHANGE_DIRECTORY;

   const std::string VALUE_MPI;
   const std::string VALUE_MPI_SINGLE;
   const std::string VALUE_SOCKETS;

   const std::string VALUE_GATHER_SCATTER;
   const std::string VALUE_POINT_TO_POINT;

   std::vector<M2NTuple> _m2ns;

   void checkDuplicates (
     const std::string& from,
     const std::string& to );
};

}} // namespace precice, m2n
