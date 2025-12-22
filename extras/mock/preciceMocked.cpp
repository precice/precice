#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <random>
#include <string>

#include <libxml/SAX2.h>

#include <precice/Exceptions.hpp>
#include <precice/Tooling.hpp>
#include <precice/Types.hpp>
#include <precice/span.hpp>
#include "utils/fmt.hpp"
// project dependencies, used to remove reliance an std types maybe not available on older versions of cpp or outside libs

// --include precise errors instead of assertions !! check when error/exception should be used
// test for exception handling (exceptions/errors when meant to be thrown)
// --mimick implicit coupling with checkpointing
// catch errors from course !! do course against mock
// replace python (and fenics) adapter mocks
// check if vector components are written correctly (later/put in issues)
// --read config (parts) to see if adapter does what config expects from it !!
// --names, dims...
// --libxml2 wrapper in utils or xml !! Not feasible due to requiring linking the full precice library
// allow optional file input for mock config path do data for now no data reading !!(or modify precice config)
// data from file exists as ASTE already

// --figure out why compilation was weird (weird issue when creating the current git structure merged old version of makefile incorrectly)
// --Ltrans warning is not easily avoidable (has to do with the way linking is handled in the entire project)
// --only happens on build from a clean state

// Add dummy data exchange
// user specifies timeseries of data from file
// run an open foam example against the mock (flap has this data alt. flow over heated plate)

// dont use random data instead use linear transform

// --add pre-commit

// ssh link for git

// add documentation to classes and functions

// add example data to mock directory and maybe integratiuon test

// --use assertions where applicable

// API declaration
namespace precice {

// Using definition from Participant.hpp
// IMPORTANT: Must match the real preCICE definition for ABI compatibility
// when precice moves to C++20, this can be replaced with std::string_view and std::span and include <string_view> and <span>
using string_view = ::precice::span<const char>;

// Forward declaration of implementation class to hide implementation details and remove need for a header
namespace impl {
class ParticipantImpl;
}

class Participant {
public:
  Participant(::precice::string_view participantName,
              ::precice::string_view configurationFileName,
              int                    solverProcessIndex,
              int                    solverProcessSize);

  Participant(::precice::string_view participantName,
              ::precice::string_view configurationFileName,
              int                    solverProcessIndex,
              int                    solverProcessSize,
              void                  *communicator);

  ~Participant();

  // Steering methods
  void initialize();
  void advance(double computedTimeStepSize);
  void finalize();

  // Implicit coupling
  bool requiresWritingCheckpoint();
  bool requiresReadingCheckpoint();

  // Status queries
  int    getMeshDimensions(::precice::string_view meshName) const;
  int    getDataDimensions(::precice::string_view meshName,
                           ::precice::string_view dataName) const;
  bool   isCouplingOngoing() const;
  bool   isTimeWindowComplete() const;
  double getMaxTimeStepSize() const;

  // Mesh access
  bool     requiresMeshConnectivityFor(::precice::string_view meshName) const;
  void     resetMesh(::precice::string_view meshName);
  VertexID setMeshVertex(::precice::string_view      meshName,
                         precice::span<const double> position);
  int      getMeshVertexSize(::precice::string_view meshName) const;
  void     setMeshVertices(::precice::string_view      meshName,
                           precice::span<const double> coordinates,
                           precice::span<VertexID>     ids);
  void     setMeshEdge(::precice::string_view meshName, VertexID first, VertexID second);
  void     setMeshEdges(::precice::string_view        meshName,
                        precice::span<const VertexID> ids);
  void     setMeshTriangle(::precice::string_view meshName, VertexID first, VertexID second, VertexID third);
  void     setMeshTriangles(::precice::string_view        meshName,
                            precice::span<const VertexID> ids);
  void     setMeshQuad(::precice::string_view meshName, VertexID first, VertexID second, VertexID third, VertexID fourth);
  void     setMeshQuads(::precice::string_view        meshName,
                        precice::span<const VertexID> ids);
  void     setMeshTetrahedron(::precice::string_view meshName, VertexID first, VertexID second, VertexID third, VertexID fourth);
  void     setMeshTetrahedra(::precice::string_view        meshName,
                             precice::span<const VertexID> ids);

  // Data access
  bool requiresInitialData();
  void writeData(::precice::string_view        meshName,
                 ::precice::string_view        dataName,
                 precice::span<const VertexID> ids,
                 precice::span<const double>   values);
  void readData(::precice::string_view        meshName,
                ::precice::string_view        dataName,
                precice::span<const VertexID> ids,
                double                        relativeReadTime,
                precice::span<double>         values) const;

  // Just-in-time mapping (experimental)
  void writeAndMapData(::precice::string_view      meshName,
                       ::precice::string_view      dataName,
                       precice::span<const double> coordinates,
                       precice::span<const double> values);
  void mapAndReadData(::precice::string_view      meshName,
                      ::precice::string_view      dataName,
                      precice::span<const double> coordinates,
                      double                      relativeReadTime,
                      precice::span<double>       values) const;

  // Direct access
  void setMeshAccessRegion(::precice::string_view      meshName,
                           precice::span<const double> boundingBox) const;
  void getMeshVertexIDsAndCoordinates(::precice::string_view  meshName,
                                      precice::span<VertexID> ids,
                                      precice::span<double>   coordinates) const;

  // Gradient data (experimental)
  bool requiresGradientDataFor(::precice::string_view meshName,
                               ::precice::string_view dataName) const;
  void writeGradientData(::precice::string_view        meshName,
                         ::precice::string_view        dataName,
                         precice::span<const VertexID> ids,
                         precice::span<const double>   gradients);

  // Profiling
  void startProfilingSection(::precice::string_view sectionName);
  void stopLastProfilingSection();

private:
  std::unique_ptr<impl::ParticipantImpl> _impl;
};

} // namespace precice

/*Implementation
  This is the implementation of the precice mock Participant class.
  It provides a minimal mock implementation of the preCICE Participant API,
  suitable for testing and demonstration purposes.
  Currently only basic checks are performed.
  Currently, data reading returns seeded generate random data,
  using data from a user-specified file if provided at initialization is currently being implemented.
  */
namespace precice {

namespace impl {

class ParticipantImpl {
public:
  ParticipantImpl(::precice::string_view participantName,
                  ::precice::string_view configurationFileName,
                  int                    solverProcessIndex,
                  int                    solverProcessSize,
                  void                  *communicator = nullptr);

  ~ParticipantImpl() = default;

  // Configuration data structures
  struct MeshInfo {
    std::string name;
    int         dimensions = 3;
    bool        provided   = false;
    bool        received   = false;
  };

  struct DataInfo {
    std::string name;
    std::string meshName;
    int         dimensions = 1;
    bool        isRead     = false;
    bool        isWrite    = false;
  };

  struct ConfigData {
    std::map<std::string, MeshInfo> meshes;
    std::vector<DataInfo>           dataItems;
    bool                            isImplicitCoupling    = false;
    bool                            hasConvergenceMeasure = false;
    int                             maxIterations         = 1;
    int                             currentIteration      = 0;
    bool                            iterationConverged    = false;
  };

  enum class DataMode {
    Random,      // Mode 1: Random data
    Buffer,      // Mode 2: Return written data buffer (default)
    ScaledBuffer // Mode 3: Return written data buffer scaled/multiplied
  };

  struct MockDataConfig {
    std::string         dataName;
    std::string         meshName;
    DataMode            mode             = DataMode::Buffer; // Default to buffer mode
    double              scalarMultiplier = 1.0;
    std::vector<double> vectorMultiplier; // Empty means use scalar
  };

  struct MockConfig {
    std::map<std::string, MockDataConfig> dataConfigs; // Key: "meshName:dataName"
    DataMode                              defaultMode             = DataMode::Buffer;
    double                                defaultScalarMultiplier = 1.0;
    std::vector<double>                   defaultVectorMultiplier; // Empty means use scalar
  };

  // Minimal internal state
  std::string                  name;
  std::string                  config;
  std::string                  mockConfigPath;
  int                          rank            = 0;
  bool                         primary         = false;
  int                          size            = 1;
  void                        *comm            = nullptr;
  bool                         initialized     = false;
  bool                         finalized       = false;
  bool                         couplingOngoing = false;
  double                       maxTimeStep     = 0.0;
  mutable std::recursive_mutex mtx;

  uint32_t    seed        = 0;
  std::size_t currentStep = 0;
  double      currentTime = 0.0;

  ConfigData configData;
  bool       configParsed = false;

  MockConfig mockConfig;
  bool       mockConfigParsed = false;
  bool       mockConfigExists = false;

  // Write data buffers: Key is "meshName:dataName"
  std::map<std::string, std::vector<double>> writeBuffers;
  // Vertex counts per mesh to detect empty provided meshes
  std::map<std::string, std::size_t> meshVertexCounts;

  enum class DataType { Scalar = 1,
                        Vector = 3 };

  // SAX parsing state for config
  struct ConfigParseState {
    bool                            inParticipant = false;
    std::string                     participantName;
    std::map<std::string, DataType> dataType; // type marker: scalar or vector
  } configParseState;

  // SAX parsing state for mock config
  struct MockParseState {
    bool                inMockedData = false;
    bool                inDefault    = false;
    std::string         currentMesh;
    std::string         currentData;
    DataMode            currentMode   = DataMode::Buffer;
    double              currentScalar = 1.0;
    std::vector<double> currentVector;
  } mockParseState;

  // Configuration parsing methods
  void parseConfig();
  void parseMockConfig();

  // Helper method for XML parsing
  void parseXMLFile(const std::string &filePath,
                    void (*startHandler)(void *, const xmlChar *, const xmlChar *, const xmlChar *, int, const xmlChar **, int, int, const xmlChar **),
                    void (*endHandler)(void *, const xmlChar *, const xmlChar *, const xmlChar *));

  // SAX callback handlers (static functions for libxml2)
  static void onConfigStartElement(void *ctx, const xmlChar *localname, const xmlChar *prefix,
                                   const xmlChar *URI, int nb_namespaces, const xmlChar **namespaces,
                                   int nb_attributes, int nb_defaulted, const xmlChar **attributes);
  static void onConfigEndElement(void *ctx, const xmlChar *localname, const xmlChar *prefix, const xmlChar *URI);
  static void onMockStartElement(void *ctx, const xmlChar *localname, const xmlChar *prefix,
                                 const xmlChar *URI, int nb_namespaces, const xmlChar **namespaces,
                                 int nb_attributes, int nb_defaulted, const xmlChar **attributes);
  static void onMockEndElement(void *ctx, const xmlChar *localname, const xmlChar *prefix, const xmlChar *URI);
};

} // namespace impl

// ParticipantImpl constructor definition
impl::ParticipantImpl::ParticipantImpl(::precice::string_view participantName,
                                       ::precice::string_view configurationFileName,
                                       int                    solverProcessIndex,
                                       int                    solverProcessSize,
                                       void                  *communicator)
    : name(std::string(participantName.data(), participantName.size())),
      config(std::string(configurationFileName.data(), configurationFileName.size())),
      rank(solverProcessIndex),
      primary(rank == 0),
      size(solverProcessSize),
      comm(communicator)
{
  // Basic invariant checks for the mock participant
  if (participantName.empty()) {
    throw precice::Error("Participant name is empty.");
  }
  if (solverProcessSize <= 0) {
    throw precice::Error(precice::utils::format_or_error("solverProcessSize must be > 0, but is {}.", solverProcessSize));
  }
  if (solverProcessIndex < 0) {
    throw precice::Error(precice::utils::format_or_error("solverProcessIndex must be >= 0, but is {}.", solverProcessIndex));
  }
  if (solverProcessIndex >= solverProcessSize) {
    throw precice::Error(precice::utils::format_or_error("solverProcessIndex={} must be smaller than solverProcessSize={}.", solverProcessIndex, solverProcessSize));
  }
  seed = static_cast<uint32_t>(rank) ^ 0x9e3779b9u;
}

// SAX callback handlers for config parsing
void impl::ParticipantImpl::onConfigStartElement(void *ctx, const xmlChar *localname, const xmlChar *prefix,
                                                 const xmlChar *URI, int nb_namespaces, const xmlChar **namespaces,
                                                 int nb_attributes, int nb_defaulted, const xmlChar **attributes)
{
  auto       *impl = static_cast<impl::ParticipantImpl *>(ctx);
  std::string elemName(reinterpret_cast<const char *>(localname));
  std::string nsPrefix = prefix ? std::string(reinterpret_cast<const char *>(prefix)) : "";

  // Parse attributes into map
  std::map<std::string, std::string> attrs;
  unsigned int                       index = 0;
  for (int i = 0; i < nb_attributes; ++i, index += 5) {
    std::string attrName(reinterpret_cast<const char *>(attributes[index]));
    auto        valueBegin = reinterpret_cast<const char *>(attributes[index + 3]);
    auto        valueEnd   = reinterpret_cast<const char *>(attributes[index + 4]);
    attrs[attrName]        = std::string(valueBegin, valueEnd - valueBegin);
  }

  auto getAttr = [&attrs](const std::string &name) -> std::string {
    auto it = attrs.find(name);
    return it != attrs.end() ? it->second : "";
  };

  auto attrIsTrue = [&attrs](const std::string &name) -> bool {
    auto it = attrs.find(name);
    if (it == attrs.end())
      return false;
    return it->second == "yes" || it->second == "true" || it->second == "1";
  };

  // Handle data type declarations
  if (nsPrefix == "data" && elemName == "scalar") {
    std::string dataName = getAttr("name");
    if (!dataName.empty())
      impl->configParseState.dataType[dataName] = DataType::Scalar;
  } else if (nsPrefix == "data" && elemName == "vector") {
    std::string dataName = getAttr("name");
    if (!dataName.empty())
      impl->configParseState.dataType[dataName] = DataType::Vector;
  }
  // Handle mesh declarations
  else if (elemName == "mesh") {
    std::string meshName = getAttr("name");
    if (!meshName.empty()) {
      impl::ParticipantImpl::MeshInfo meshInfo;
      meshInfo.name                     = meshName;
      std::string dimsStr               = getAttr("dimensions");
      meshInfo.dimensions               = dimsStr.empty() ? 3 : std::stoi(dimsStr);
      meshInfo.provided                 = attrIsTrue("provide");
      meshInfo.received                 = attrIsTrue("receive");
      impl->configData.meshes[meshName] = meshInfo;
    }
  }
  // Handle participant
  else if (elemName == "participant") {
    std::string participantName = getAttr("name");
    if (participantName == impl->name) {
      impl->configParseState.inParticipant   = true;
      impl->configParseState.participantName = participantName;
    }
  }
  // Inside participant
  else if (impl->configParseState.inParticipant) {
    if (elemName == "provide-mesh" || elemName == "receive-mesh" || elemName == "use-mesh") {
      std::string meshName = getAttr("name");
      auto        it       = impl->configData.meshes.find(meshName);
      if (it != impl->configData.meshes.end()) {
        if (elemName == "provide-mesh" || attrIsTrue("provide"))
          it->second.provided = true;
        if (elemName == "receive-mesh" || attrIsTrue("receive"))
          it->second.received = true;
      }
    } else if (elemName == "write-data") {
      std::string dataName = getAttr("name");
      std::string meshName = getAttr("mesh");
      if (!dataName.empty() && !meshName.empty()) {
        impl::ParticipantImpl::DataInfo dataInfo;
        dataInfo.name     = dataName;
        dataInfo.meshName = meshName;
        dataInfo.isWrite  = true;
        // Determine data dimensions: for vectors, use mesh dimensions; for scalars, use 1
        auto typeIt   = impl->configParseState.dataType.find(dataName);
        bool isVector = (typeIt != impl->configParseState.dataType.end() && typeIt->second == DataType::Vector);
        if (isVector) {
          // Vector data: use mesh dimensions
          auto meshIt         = impl->configData.meshes.find(meshName);
          dataInfo.dimensions = (meshIt != impl->configData.meshes.end()) ? meshIt->second.dimensions : 3;
        } else {
          // Scalar data: always 1D
          dataInfo.dimensions = 1;
        }
        impl->configData.dataItems.push_back(dataInfo);
      }
    } else if (elemName == "read-data") {
      std::string dataName = getAttr("name");
      std::string meshName = getAttr("mesh");
      if (!dataName.empty() && !meshName.empty()) {
        impl::ParticipantImpl::DataInfo dataInfo;
        dataInfo.name     = dataName;
        dataInfo.meshName = meshName;
        dataInfo.isRead   = true;
        // Determine data dimensions: for vectors, use mesh dimensions; for scalars, use 1
        auto typeIt   = impl->configParseState.dataType.find(dataName);
        bool isVector = (typeIt != impl->configParseState.dataType.end() && typeIt->second == DataType::Vector);
        if (isVector) {
          // Vector data: use mesh dimensions
          auto meshIt         = impl->configData.meshes.find(meshName);
          dataInfo.dimensions = (meshIt != impl->configData.meshes.end()) ? meshIt->second.dimensions : 3;
        } else {
          // Scalar data: always 1D
          dataInfo.dimensions = 1;
        }
        impl->configData.dataItems.push_back(dataInfo);
      }
    }
  }
  // Handle coupling schemes
  if (nsPrefix == "coupling-scheme" && (elemName == "serial-implicit" || elemName == "parallel-implicit")) {
    impl->configData.isImplicitCoupling = true;
  } else if (impl->configData.isImplicitCoupling && elemName == "max-iterations") {
    std::string maxIterStr = getAttr("value");
    if (!maxIterStr.empty()) {
      impl->configData.maxIterations = std::stoi(maxIterStr);
    }
  } else if (impl->configData.isImplicitCoupling &&
             (elemName == "relative-convergence-measure" || elemName == "absolute-convergence-measure" ||
              elemName == "residual-relative-convergence-measure" || elemName == "min-iteration-convergence-measure")) {
    impl->configData.hasConvergenceMeasure = true;
  }
}

void impl::ParticipantImpl::onConfigEndElement(void *ctx, const xmlChar *localname, const xmlChar *prefix, const xmlChar *URI)
{
  auto       *impl = static_cast<impl::ParticipantImpl *>(ctx);
  std::string elemName(reinterpret_cast<const char *>(localname));
  if (elemName == "participant") {
    impl->configParseState.inParticipant = false;
  }
}

// Helper method to parse XML files using SAX
void impl::ParticipantImpl::parseXMLFile(const std::string &filePath,
                                         void (*startHandler)(void *, const xmlChar *, const xmlChar *, const xmlChar *, int, const xmlChar **, int, int, const xmlChar **),
                                         void (*endHandler)(void *, const xmlChar *, const xmlChar *, const xmlChar *))
{
  // Setup SAX handler
  xmlSAXHandler saxHandler;
  std::memset(&saxHandler, 0, sizeof(xmlSAXHandler));
  saxHandler.initialized    = XML_SAX2_MAGIC;
  saxHandler.startElementNs = startHandler;
  saxHandler.endElementNs   = endHandler;

  // Read file content
  std::ifstream ifs(filePath);
  if (!ifs) {
    throw precice::Error(precice::utils::format_or_error(
        "Could not open configuration file '{}'", filePath));
  }
  std::string content{std::istreambuf_iterator<char>(ifs), std::istreambuf_iterator<char>()};
  if (content.empty()) {
    throw precice::Error(precice::utils::format_or_error(
        "Configuration file '{}' is empty", filePath));
  }

  // Parse using SAX
  xmlParserCtxtPtr ctxt = xmlCreatePushParserCtxt(&saxHandler, static_cast<void *>(this),
                                                  content.c_str(), content.size(), nullptr);
  xmlParseChunk(ctxt, nullptr, 0, 1);
  xmlFreeParserCtxt(ctxt);
  xmlCleanupParser();
}

void impl::ParticipantImpl::parseConfig()
{
  if (configParsed)
    return;

  try {
    // Parse configuration file
    parseXMLFile(config, onConfigStartElement, onConfigEndElement);

    if (!configParseState.inParticipant && configParseState.participantName.empty()) {
      throw precice::Error(precice::utils::format_or_error(
          "Participant '{}' not found in configuration '{}'", name, config));
    }

    // If no mesh was explicitly marked as provided but meshes exist, assume all non-received meshes are provided (fallback for simplified configs)
    bool anyProvided = false;
    for (const auto &m : configData.meshes) {
      if (m.second.provided) {
        anyProvided = true;
        break;
      }
    }
    if (!anyProvided) {
      for (auto &m : configData.meshes) {
        if (!m.second.received) {
          m.second.provided = true;
        }
      }
    }

    // Ensure meshes referenced by data items exist
    for (const auto &dataInfo : configData.dataItems) {
      auto it = configData.meshes.find(dataInfo.meshName);
      if (it == configData.meshes.end()) {
        impl::ParticipantImpl::MeshInfo mi;
        mi.name                              = dataInfo.meshName;
        mi.dimensions                        = 3;
        mi.provided                          = dataInfo.isWrite;
        mi.received                          = dataInfo.isRead;
        configData.meshes[dataInfo.meshName] = mi;
      } else {
        if (dataInfo.isWrite)
          it->second.provided = true;
        if (dataInfo.isRead)
          it->second.received = true;
      }
    }

    // Set defaults for implicit coupling
    if (configData.isImplicitCoupling) {
      // Only print warning if no convergence measure was found
      if (!configData.hasConvergenceMeasure) {
        std::cerr << "---[precice]  WARNING: No convergence measures were defined for an implicit coupling scheme. "
                  << "It will always iterate the maximum amount iterations, which is " << configData.maxIterations
                  << ". You may want to add a convergence measure in your <coupling-scheme:.../> in your configuration."
                  << std::endl;
      }
    }

    configParsed = true;
  } catch (const precice::Error &) {
    throw;
  } catch (const std::exception &e) {
    std::cerr << "Warning: Could not parse configuration file '" << config
              << "': " << e.what() << ". Using limited mock behavior." << std::endl;
    configParsed = true;
  }
}

// SAX callback handlers for mock config parsing
void impl::ParticipantImpl::onMockStartElement(void *ctx, const xmlChar *localname, const xmlChar *prefix,
                                               const xmlChar *URI, int nb_namespaces, const xmlChar **namespaces,
                                               int nb_attributes, int nb_defaulted, const xmlChar **attributes)
{
  auto       *impl = static_cast<impl::ParticipantImpl *>(ctx);
  std::string elemName(reinterpret_cast<const char *>(localname));

  // Parse attributes into map
  std::map<std::string, std::string> attrs;
  unsigned int                       index = 0;
  for (int i = 0; i < nb_attributes; ++i, index += 5) {
    std::string attrName(reinterpret_cast<const char *>(attributes[index]));
    auto        valueBegin = reinterpret_cast<const char *>(attributes[index + 3]);
    auto        valueEnd   = reinterpret_cast<const char *>(attributes[index + 4]);
    attrs[attrName]        = std::string(valueBegin, valueEnd - valueBegin);
  }

  if (elemName == "default") {
    impl->mockParseState.inDefault = true;
    std::string modeStr            = attrs["mode"];
    if (modeStr == "random") {
      impl->mockConfig.defaultMode = DataMode::Random;
    } else if (modeStr == "scaled") {
      impl->mockConfig.defaultMode = DataMode::ScaledBuffer;
    } else {
      impl->mockConfig.defaultMode = DataMode::Buffer;
    }
    impl->mockParseState.currentScalar = 1.0;
    impl->mockParseState.currentVector.clear();
  } else if (elemName == "mocked-data") {
    impl->mockParseState.inMockedData = true;
    impl->mockParseState.currentMesh  = attrs["mesh"];
    impl->mockParseState.currentData  = attrs["data"];
    std::string modeStr               = attrs["mode"];
    if (modeStr == "random") {
      impl->mockParseState.currentMode = DataMode::Random;
    } else if (modeStr == "scaled") {
      impl->mockParseState.currentMode = DataMode::ScaledBuffer;
    } else {
      impl->mockParseState.currentMode = DataMode::Buffer;
    }
    impl->mockParseState.currentScalar = 1.0;
    impl->mockParseState.currentVector.clear();
  } else if (impl->mockParseState.inMockedData || impl->mockParseState.inDefault) {
    if (elemName == "scalar-multiplier") {
      std::string valueStr = attrs["value"];
      if (!valueStr.empty()) {
        impl->mockParseState.currentScalar = std::stod(valueStr);
      }
    } else if (elemName == "vector-multiplier") {
      std::string valuesStr = attrs["values"];
      if (!valuesStr.empty()) {
        size_t pos = 0;
        while (pos < valuesStr.length()) {
          size_t nextSemi = valuesStr.find(';', pos);
          if (nextSemi == std::string::npos)
            nextSemi = valuesStr.length();
          std::string valStr = valuesStr.substr(pos, nextSemi - pos);
          if (!valStr.empty())
            impl->mockParseState.currentVector.push_back(std::stod(valStr));
          pos = nextSemi + 1;
        }
      }
    }
  }
}

void impl::ParticipantImpl::onMockEndElement(void *ctx, const xmlChar *localname, const xmlChar *prefix, const xmlChar *URI)
{
  auto       *impl = static_cast<impl::ParticipantImpl *>(ctx);
  std::string elemName(reinterpret_cast<const char *>(localname));
  if (elemName == "default" && impl->mockParseState.inDefault) {
    // Apply default multipliers
    impl->mockConfig.defaultScalarMultiplier = impl->mockParseState.currentScalar;
    impl->mockConfig.defaultVectorMultiplier = impl->mockParseState.currentVector;
    impl->mockParseState.inDefault           = false;
  } else if (elemName == "mocked-data" && impl->mockParseState.inMockedData) {
    if (!impl->mockParseState.currentMesh.empty() && !impl->mockParseState.currentData.empty()) {
      std::string    key = impl->mockParseState.currentMesh + ":" + impl->mockParseState.currentData;
      MockDataConfig config;
      config.meshName                   = impl->mockParseState.currentMesh;
      config.dataName                   = impl->mockParseState.currentData;
      config.mode                       = impl->mockParseState.currentMode;
      config.scalarMultiplier           = impl->mockParseState.currentScalar;
      config.vectorMultiplier           = impl->mockParseState.currentVector;
      impl->mockConfig.dataConfigs[key] = config;
    }
    impl->mockParseState.inMockedData = false;
  }
}

void impl::ParticipantImpl::parseMockConfig()
{
  if (mockConfigParsed)
    return;

  try {
    // Derive mock-config path from precice config path
    mockConfigPath   = config;
    size_t lastSlash = mockConfigPath.find_last_of("/\\");
    size_t lastDot   = mockConfigPath.find_last_of('.');

    if (lastSlash != std::string::npos) {
      if (lastDot != std::string::npos && lastDot > lastSlash) {
        mockConfigPath = mockConfigPath.substr(0, lastSlash + 1) + "mock-config.xml";
      } else {
        mockConfigPath = mockConfigPath + "/mock-config.xml";
      }
    } else {
      // No directory, same directory as executable
      mockConfigPath = "mock-config.xml";
    }

    // Check if mock config file exists
    std::ifstream mockFile(mockConfigPath);
    if (!mockFile.good()) {
      // No mock config - use defaults
      mockConfig.defaultMode = DataMode::Buffer;
      mockConfigParsed       = true;
      mockConfigExists       = false;
      return;
    }
    mockFile.close();
    mockConfigExists = true;

    // Parse mock configuration file
    parseXMLFile(mockConfigPath, onMockStartElement, onMockEndElement);
    mockConfigParsed = true;

  } catch (const std::exception &e) {
    std::cerr << "Warning: Could not parse mock configuration file: " << e.what()
              << ". Using default buffer mode." << std::endl;
    mockConfig.defaultMode = DataMode::Buffer;
    mockConfigParsed       = true;
  }
}

// Participant constructors / destructor

Participant::Participant(
    ::precice::string_view participantName,
    ::precice::string_view configurationFileName,
    int                    solverProcessIndex,
    int                    solverProcessSize)
    : _impl(std::make_unique<impl::ParticipantImpl>(participantName,
                                                    configurationFileName,
                                                    solverProcessIndex,
                                                    solverProcessSize,
                                                    nullptr))
{
  // Parse configuration at construction time so that mesh/data queries
  // (and bindings) can work before initialize() is called — matching
  // behavior of the real preCICE library.
  _impl->parseConfig();
  _impl->parseMockConfig();
}

Participant::Participant(
    ::precice::string_view participantName,
    ::precice::string_view configurationFileName,
    int                    solverProcessIndex,
    int                    solverProcessSize,
    void                  *communicator)
    : _impl(std::make_unique<impl::ParticipantImpl>(participantName,
                                                    configurationFileName,
                                                    solverProcessIndex,
                                                    solverProcessSize,
                                                    communicator))
{
  // Parse configuration at construction time so that mesh/data queries
  // (and bindings) can work before initialize() is called — matching
  // behavior of the real preCICE library.
  _impl->parseConfig();
  _impl->parseMockConfig();
}

Participant::~Participant() = default;

// Steering methods

void Participant::initialize()
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (_impl->finalized) {
    throw precice::Error("initialize() cannot be called after finalize().");
  }
  if (_impl->initialized) {
    throw precice::Error("initialize() may only be called once.");
  }

  // Print startup banner (only on rank 0)
  if (_impl->primary) {
    std::cout << "---[precice]  This is preCICE version 3.3.0 (mock)" << std::endl;
    std::cout << "---[precice]  Revision info: mock-implementation" << std::endl;
    std::cout << "---[precice]  Build type: Release (mock)" << std::endl;
    try {
      std::cout << "---[precice]  Working directory \"" << std::filesystem::current_path().string() << "\"" << std::endl;
    } catch (std::filesystem::filesystem_error &fse) {
      std::cout << "---[precice]  Working directory unknown due to error \"" << fse.what() << "\"" << std::endl;
    }
    if (_impl->mockConfigExists) {
      std::cout << "---[precice]  Configuring preCICE (mock) with configuration \"" << std::filesystem::absolute(_impl->config).string() << "\" and mock configuration \"" << std::filesystem::absolute(_impl->mockConfigPath).string() << "\"" << std::endl;
    } else {
      std::cout << "---[precice]  Configuring preCICE (mock) with configuration \"" << std::filesystem::absolute(_impl->config).string() << "\" and default mock configuration" << std::endl;
    }
    std::cout << "---[precice]  I am participant \"" << _impl->name << "\"" << std::endl;
  }

  // Communication setup prints
  if (_impl->primary) {
    std::cout << "---[precice]  Setting up primary communication to coupling partner/s" << std::endl;
    std::cout << "---[precice]  Primary ranks are connected" << std::endl;
    std::cout << "---[precice]  Setting up preliminary secondary communication to coupling partner/s" << std::endl;
  }

  // Mesh preparation prints and empty-mesh detection for provided meshes
  for (const auto &meshEntry : _impl->configData.meshes) {
    const auto &meshInfo = meshEntry.second;
    if (!meshInfo.provided) {
      continue;
    }

    if (_impl->primary) {
      std::cout << "---[precice]  Prepare partition for mesh " << meshInfo.name << std::endl;
      std::cout << "---[precice]  Gather mesh " << meshInfo.name << std::endl;
      std::cout << "---[precice]  Send global mesh " << meshInfo.name << std::endl;
    }

    std::size_t count = 0;
    auto        it    = _impl->meshVertexCounts.find(meshInfo.name);
    if (it != _impl->meshVertexCounts.end()) {
      count = it->second;
    }
    if (count == 0) {
      throw precice::Error(precice::utils::format_or_error(
          "The provided mesh \"{}\" is empty. Please set the mesh using setMeshVertex()/setMeshVertices() prior to calling initialize().",
          meshInfo.name));
    }
  }

  _impl->initialized     = true;
  _impl->couplingOngoing = true;
  _impl->maxTimeStep     = 1.0;
  _impl->currentStep     = 0;
  _impl->currentTime     = 0.0;

  // Initialize implicit coupling state
  if (_impl->configData.isImplicitCoupling) {
    _impl->configData.currentIteration   = 0;
    _impl->configData.iterationConverged = false;
  }
}

void Participant::advance(double computedTimeStepSize)
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (!_impl->initialized) {
    throw precice::Error("initialize() has to be called before advance().");
  }
  if (_impl->finalized) {
    throw precice::Error("advance() cannot be called after finalize().");
  }
  if (!isCouplingOngoing()) {
    throw precice::Error("advance() cannot be called when isCouplingOngoing() returns false.");
  }
  if (!std::isfinite(computedTimeStepSize)) {
    throw precice::Error("advance() cannot be called with an infinite time step size.");
  }
  if (computedTimeStepSize == 0.0) {
    throw precice::Error("advance() cannot be called with a time step size of 0.");
  }
  if (computedTimeStepSize < 0.0) {
    throw precice::Error(precice::utils::format_or_error("advance() cannot be called with a negative time step size {}.", computedTimeStepSize));
  }
  if (_impl->maxTimeStep == _impl->currentStep) {
    _impl->couplingOngoing = false;
    return;
  }

  // Handle implicit coupling iterations
  if (_impl->configData.isImplicitCoupling) {
    _impl->configData.currentIteration++;

    // Simple convergence logic: converge after max iterations
    if (_impl->configData.currentIteration >= _impl->configData.maxIterations) {
      _impl->configData.iterationConverged = true;
      _impl->configData.currentIteration   = 0;
      _impl->currentStep += 1;
      _impl->currentTime += computedTimeStepSize;
    }
  } else {
    // Explicit coupling: always advance
    _impl->currentStep += 1;
    _impl->currentTime += computedTimeStepSize;
  }

  if (computedTimeStepSize > 0.0) {
    _impl->maxTimeStep = std::max(_impl->maxTimeStep, computedTimeStepSize);
  }
}

void Participant::finalize()
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (_impl->finalized) {
    throw precice::Error("finalize() may only be called once.");
  }

  // Print finalize message (only on rank 0)
  if (_impl->primary) {
    std::cout << "---[precice]  Finalizing preCICE (mock)" << std::endl;
  }

  _impl->initialized     = false;
  _impl->couplingOngoing = false;
  _impl->finalized       = true;
}

// Implicit coupling

bool Participant::requiresWritingCheckpoint()
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (!_impl->initialized) {
    throw precice::Error("initialize() has to be called before requiresWritingCheckpoint().");
  }

  // For implicit coupling, require checkpoint at start of iteration
  if (_impl->configData.isImplicitCoupling) {
    // Checkpoint needed at beginning of time window (iteration 0)
    return _impl->configData.currentIteration == 0 && !_impl->configData.iterationConverged;
  }

  return false;
}

bool Participant::requiresReadingCheckpoint()
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (!_impl->initialized) {
    throw precice::Error("initialize() has to be called before requiresReadingCheckpoint().");
  }

  // For implicit coupling, require reading checkpoint when not converged
  if (_impl->configData.isImplicitCoupling) {
    // Read checkpoint needed when iteration hasn't converged and we're past iteration 0
    return _impl->configData.currentIteration > 0 && !_impl->configData.iterationConverged;
  }

  return false;
}

// Status queries

int Participant::getMeshDimensions(::precice::string_view meshName) const
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);

  std::string meshNameStr(meshName.data(), meshName.size());

  // Check if config was parsed and mesh exists
  if (_impl->configParsed) {
    auto it = _impl->configData.meshes.find(meshNameStr);
    if (it != _impl->configData.meshes.end()) {
      return it->second.dimensions;
    } else {
      throw precice::Error(precice::utils::format_or_error(
          "Mesh '{}' is not used by participant '{}'. Available meshes: {}",
          meshNameStr, _impl->name, [&]() {
            std::string meshes;
            for (const auto &m : _impl->configData.meshes) {
              if (!meshes.empty())
                meshes += ", ";
              meshes += m.first;
            }
            return meshes.empty() ? "none" : meshes;
          }()));
    }
  }

  // If configuration wasn't parsed or mesh wasn't found in parsed config,
  // try to return a previously inferred mesh dimension (e.g. from earlier
  // setMeshVertex/setMeshVertices calls). If none is available, return 0
  // to indicate unknown dimensionality so higher-level bindings can infer
  // dimensions from the provided vertex data (matching real preCICE behavior).
  auto it2 = _impl->configData.meshes.find(meshNameStr);
  if (it2 != _impl->configData.meshes.end()) {
    return it2->second.dimensions;
  }

  // Unknown mesh dimensionality
  return 0;
}

int Participant::getDataDimensions(::precice::string_view meshName,
                                   ::precice::string_view dataName) const
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);

  std::string meshNameStr(meshName.data(), meshName.size());
  std::string dataNameStr(dataName.data(), dataName.size());

  // Check if config was parsed and data exists
  if (_impl->configParsed) {
    for (const auto &dataInfo : _impl->configData.dataItems) {
      if (dataInfo.name == dataNameStr && dataInfo.meshName == meshNameStr) {
        return dataInfo.dimensions;
      }
    }

    // Data not found - throw error
    throw precice::Error(precice::utils::format_or_error(
        "Data '{}' on mesh '{}' is not used by participant '{}'",
        dataNameStr, meshNameStr, _impl->name));
  }

  // Fallback: mock as scalar
  return 1;
}

bool Participant::isCouplingOngoing() const
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (_impl->finalized) {
    throw precice::Error("isCouplingOngoing() cannot be called after finalize().");
  }
  if (!_impl->initialized) {
    throw precice::Error("initialize() has to be called before isCouplingOngoing() can be evaluated.");
  }
  return _impl->couplingOngoing;
}

bool Participant::isTimeWindowComplete() const
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (!_impl->initialized) {
    throw precice::Error("initialize() has to be called before isTimeWindowComplete().");
  }
  if (_impl->finalized) {
    throw precice::Error("isTimeWindowComplete() cannot be called after finalize().");
  }

  // For implicit coupling, time window is complete when converged
  if (_impl->configData.isImplicitCoupling) {
    return _impl->configData.iterationConverged;
  }

  return !_impl->couplingOngoing;
}

double Participant::getMaxTimeStepSize() const
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (_impl->finalized) {
    throw precice::Error("getMaxTimeStepSize() cannot be called after finalize().");
  }
  if (!_impl->initialized) {
    throw precice::Error("initialize() has to be called before getMaxTimeStepSize() can be evaluated.");
  }
  return _impl->maxTimeStep;
}

// Mesh access

bool Participant::requiresMeshConnectivityFor(::precice::string_view /*meshName*/) const
{
  return false;
}

void Participant::resetMesh(::precice::string_view meshName)
{
  if (!_impl->initialized) {
    throw precice::Error("initialize() has to be called before resetMesh().");
  }
  std::string                           meshNameStr(meshName.data(), meshName.size());
  std::lock_guard<std::recursive_mutex> lkm(_impl->mtx);
  _impl->meshVertexCounts[meshNameStr] = 0;
}

VertexID Participant::setMeshVertex(
    ::precice::string_view      meshName,
    precice::span<const double> position)
{
  std::lock_guard<std::recursive_mutex> lkm(_impl->mtx);
  if (_impl->initialized) {
    throw precice::Error("setMeshVertex() cannot be called after initialize(). Mesh modification is only allowed before calling initialize().");
  }
  std::string meshNameStr(meshName.data(), meshName.size());
  // Determine expected dimensionality from configuration if available
  int  expectedDims = 3;
  auto meshIt       = _impl->configData.meshes.find(meshNameStr);
  if (meshIt != _impl->configData.meshes.end()) {
    expectedDims = meshIt->second.dimensions;
  }
  if (position.size() != static_cast<std::size_t>(expectedDims)) {
    throw precice::Error(precice::utils::format_or_error("setMeshVertex() was called with {} coordinates, but expects {} coordinates for this mesh.", position.size(), expectedDims));
  }
  auto &count = _impl->meshVertexCounts[meshNameStr];
  ++count;
  return VertexID(count - 1);
}

int Participant::getMeshVertexSize(::precice::string_view meshName) const
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  std::string                           meshNameStr(meshName.data(), meshName.size());
  auto                                  it = _impl->meshVertexCounts.find(meshNameStr);
  if (it != _impl->meshVertexCounts.end()) {
    return static_cast<int>(it->second);
  }
  return 0;
}

void Participant::setMeshVertices(
    ::precice::string_view      meshName,
    precice::span<const double> coordinates,
    precice::span<VertexID>     ids)
{
  std::lock_guard<std::recursive_mutex> lkm(_impl->mtx);
  if (_impl->initialized) {
    throw precice::Error("setMeshVertices() cannot be called after initialize(). Mesh modification is only allowed before calling initialize().");
  }
  std::string meshNameStr(meshName.data(), meshName.size());
  // Determine expected dimensionality from configuration if available
  int  expectedDims = 3;
  auto meshIt       = _impl->configData.meshes.find(meshNameStr);
  if (meshIt != _impl->configData.meshes.end()) {
    expectedDims = meshIt->second.dimensions;
  }
  // Expect coordinates as expectedDims * N entries for N VertexIDs
  if (coordinates.size() != ids.size() * expectedDims) {
    throw precice::Error(precice::utils::format_or_error("setMeshVertices() was called with {} vertices and {} coordinates ({}D), but needs {} coordinates ({} x {}).", ids.size(), coordinates.size(), expectedDims, ids.size() * expectedDims, ids.size(), expectedDims));
  }
  _impl->meshVertexCounts[meshNameStr] += ids.size();
}

void Participant::setMeshEdge(
    ::precice::string_view /*meshName*/,
    VertexID first,
    VertexID second)
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (_impl->initialized) {
    throw precice::Error("setMeshEdge() cannot be called after initialize(). Mesh modification is only allowed before calling initialize().");
  }
  // basic check: edge endpoints should not be identical
  if (first == second) {
    throw precice::Error(precice::utils::format_or_error("setMeshEdge() was called with vertices [{}, {}], but both vertices are equal.", first, second));
  }
  // no-op
}

void Participant::setMeshEdges(
    ::precice::string_view /*meshName*/,
    precice::span<const VertexID> ids)
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (_impl->initialized) {
    throw precice::Error("setMeshEdges() cannot be called after initialize(). Mesh modification is only allowed before calling initialize().");
  }
  // expect pairs of vertex ids
  if (ids.size() % 2 != 0) {
    throw precice::Error(precice::utils::format_or_error("setMeshEdges() was called with {} vertices, but needs an even number of vertices (2 per edge).", ids.size()));
  }
  // no-op
}

void Participant::setMeshTriangle(
    ::precice::string_view /*meshName*/,
    VertexID first,
    VertexID second,
    VertexID third)
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (_impl->initialized) {
    throw precice::Error("setMeshTriangle() cannot be called after initialize(). Mesh modification is only allowed before calling initialize().");
  }
  // triangle vertices should be distinct
  if (first == second || first == third || second == third) {
    throw precice::Error(precice::utils::format_or_error("setMeshTriangle() was called with vertices [{}, {}, {}], but some vertices are equal.", first, second, third));
  }
  // no-op
}

void Participant::setMeshTriangles(
    ::precice::string_view /*meshName*/,
    precice::span<const VertexID> ids)
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (_impl->initialized) {
    throw precice::Error("setMeshTriangles() cannot be called after initialize(). Mesh modification is only allowed before calling initialize().");
  }
  // expect triples of ids
  if (ids.size() % 3 != 0) {
    throw precice::Error(precice::utils::format_or_error("setMeshTriangles() was called with {} vertices, but needs a number of vertices divisible by 3 (3 per triangle).", ids.size()));
  }
  // no-op
}

void Participant::setMeshQuad(
    ::precice::string_view /*meshName*/,
    VertexID first,
    VertexID second,
    VertexID third,
    VertexID fourth)
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (_impl->initialized) {
    throw precice::Error("setMeshQuad() cannot be called after initialize(). Mesh modification is only allowed before calling initialize().");
  }
  // quad vertices should be distinct
  if (first == second || first == third || first == fourth || second == third || second == fourth || third == fourth) {
    throw precice::Error(precice::utils::format_or_error("setMeshQuad() was called with vertices [{}, {}, {}, {}], but some vertices are equal.", first, second, third, fourth));
  }
  // no-op
}

void Participant::setMeshQuads(
    ::precice::string_view /*meshName*/,
    precice::span<const VertexID> ids)
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (_impl->initialized) {
    throw precice::Error("setMeshQuads() cannot be called after initialize(). Mesh modification is only allowed before calling initialize().");
  }
  // expect groups of 4 vertex ids
  if (ids.size() % 4 != 0) {
    throw precice::Error(precice::utils::format_or_error("setMeshQuads() was called with {} vertices, but needs a number of vertices divisible by 4 (4 per quad).", ids.size()));
  }
  // no-op
}

void Participant::setMeshTetrahedron(
    ::precice::string_view /*meshName*/,
    VertexID first,
    VertexID second,
    VertexID third,
    VertexID fourth)
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (_impl->initialized) {
    throw precice::Error("setMeshTetrahedron() cannot be called after initialize(). Mesh modification is only allowed before calling initialize().");
  }
  // tetrahedron vertices should be distinct
  if (first == second || first == third || first == fourth || second == third || second == fourth || third == fourth) {
    throw precice::Error(precice::utils::format_or_error("setMeshTetrahedron() was called with vertices [{}, {}, {}, {}], but some vertices are equal.", first, second, third, fourth));
  }
  // no-op
}

void Participant::setMeshTetrahedra(
    ::precice::string_view /*meshName*/,
    precice::span<const VertexID> ids)
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (_impl->initialized) {
    throw precice::Error("setMeshTetrahedra() cannot be called after initialize(). Mesh modification is only allowed before calling initialize().");
  }
  // expect groups of 4 vertex ids
  if (ids.size() % 4 != 0) {
    throw precice::Error(precice::utils::format_or_error("setMeshTetrahedra() was called with {} vertices, but needs a number of vertices divisible by 4 (4 per tetrahedron).", ids.size()));
  }
  // no-op
}

// Data access

bool Participant::requiresInitialData()
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (_impl->initialized) {
    throw precice::Error("requiresInitialData() has to be called before initialize().");
  }
  return false;
}

void Participant::writeData(
    ::precice::string_view        meshName,
    ::precice::string_view        dataName,
    precice::span<const VertexID> ids,
    precice::span<const double>   values)
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (!_impl->initialized) {
    throw precice::Error("initialize() has to be called before writeData().");
  }

  std::string meshNameStr(meshName.data(), meshName.size());
  std::string dataNameStr(dataName.data(), dataName.size());

  // Validate against config if parsed
  if (_impl->configParsed) {
    bool found        = false;
    int  expectedDims = 1;

    for (const auto &dataInfo : _impl->configData.dataItems) {
      if (dataInfo.name == dataNameStr && dataInfo.meshName == meshNameStr) {
        if (!dataInfo.isWrite) {
          throw precice::Error(precice::utils::format_or_error(
              "Data '{}' on mesh '{}' is not configured for writing by participant '{}'. "
              "Please add <write-data name=\"{}\" mesh=\"{}\" /> to the configuration.",
              dataNameStr, meshNameStr, _impl->name, dataNameStr, meshNameStr));
        }
        found        = true;
        expectedDims = dataInfo.dimensions;
        break;
      }
    }

    if (!found) {
      throw precice::Error(precice::utils::format_or_error(
          "Data '{}' on mesh '{}' is not used by participant '{}'",
          dataNameStr, meshNameStr, _impl->name));
    }

    // Check value count matches data dimensions
    if (values.size() != ids.size() * expectedDims) {
      throw precice::Error(precice::utils::format_or_error(
          "writeData() was called with {} vertices and {} values for {}-dimensional data '{}', "
          "but needs {} values ({} x {}).",
          ids.size(), values.size(), expectedDims, dataNameStr,
          ids.size() * expectedDims, ids.size(), expectedDims));
    }
  } else {
    // Fallback validation without config
    if (values.size() != ids.size()) {
      throw precice::Error(precice::utils::format_or_error(
          "writeData() was called with {} vertices and {} values, but the number of vertices and values must match.",
          ids.size(), values.size()));
    }
  }

  // Store written data in buffer for later reading
  std::string key = meshNameStr + ":" + dataNameStr;
  _impl->writeBuffers[key].assign(values.begin(), values.end());
}

void Participant::readData(
    ::precice::string_view        meshName,
    ::precice::string_view        dataName,
    precice::span<const VertexID> ids,
    double /*relativeReadTime*/,
    precice::span<double> values) const
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  const std::size_t                     n = values.size();
  if (!_impl->initialized) {
    throw precice::Error("initialize() has to be called before readData().");
  }

  std::string meshNameStr(meshName.data(), meshName.size());
  std::string dataNameStr(dataName.data(), dataName.size());

  // Validate against config if parsed
  if (_impl->configParsed) {
    bool found        = false;
    int  expectedDims = 1;

    for (const auto &dataInfo : _impl->configData.dataItems) {
      if (dataInfo.name == dataNameStr && dataInfo.meshName == meshNameStr) {
        if (!dataInfo.isRead) {
          throw precice::Error(precice::utils::format_or_error(
              "Data '{}' on mesh '{}' is not configured for reading by participant '{}'. "
              "Please add <read-data name=\"{}\" mesh=\"{}\" /> to the configuration.",
              dataNameStr, meshNameStr, _impl->name, dataNameStr, meshNameStr));
        }
        found        = true;
        expectedDims = dataInfo.dimensions;
        break;
      }
    }

    if (!found) {
      throw precice::Error(precice::utils::format_or_error(
          "Data '{}' on mesh '{}' is not used by participant '{}'",
          dataNameStr, meshNameStr, _impl->name));
    }

    // Check value count matches data dimensions
    if (n != ids.size() * expectedDims) {
      throw precice::Error(precice::utils::format_or_error(
          "readData() was called with {} vertices and {} values for {}-dimensional data '{}', "
          "but needs {} values ({} x {}).",
          ids.size(), values.size(), expectedDims, dataNameStr,
          ids.size() * expectedDims, ids.size(), expectedDims));
    }
  } else {
    // Fallback validation without config
    if (n != ids.size()) {
      throw precice::Error(precice::utils::format_or_error(
          "readData() was called with {} vertices and {} values, but the number of vertices and values must match.",
          ids.size(), values.size()));
    }
  }

  // Determine data mode from mock config
  std::string                                  key            = meshNameStr + ":" + dataNameStr;
  impl::ParticipantImpl::DataMode              mode           = _impl->mockConfig.defaultMode;
  const impl::ParticipantImpl::MockDataConfig *mockDataConfig = nullptr;
  double                                       scalarMult     = _impl->mockConfig.defaultScalarMultiplier;
  std::vector<double>                          vectorMult     = _impl->mockConfig.defaultVectorMultiplier;

  auto mockIt = _impl->mockConfig.dataConfigs.find(key);
  if (mockIt != _impl->mockConfig.dataConfigs.end()) {
    mode           = mockIt->second.mode;
    mockDataConfig = &mockIt->second;
    scalarMult     = mockIt->second.scalarMultiplier;
    vectorMult     = mockIt->second.vectorMultiplier;
  }

  // Apply the selected mode
  switch (mode) {
  case impl::ParticipantImpl::DataMode::Random: {
    // Mode 1: Random data
    std::mt19937                           gen(static_cast<uint32_t>(_impl->seed + static_cast<uint32_t>(_impl->currentStep)));
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (std::size_t i = 0; i < n; ++i) {
      values[i] = dist(gen);
    }
    break;
  }

  case impl::ParticipantImpl::DataMode::Buffer: {
    // Mode 2: Return buffered write data
    auto bufferIt = _impl->writeBuffers.find(key);
    if (bufferIt != _impl->writeBuffers.end() && !bufferIt->second.empty()) {
      const auto &buffer = bufferIt->second;
      // Strict validation: buffer size must match requested size
      if (buffer.size() != n) {
        throw precice::Error(precice::utils::format_or_error(
            "readData() was called with {} values for data '{}' on mesh '{}', "
            "but writeData() previously wrote {} values. "
            "The number of values in readData() must match writeData() for the same mesh/data pair. "
            "This typically indicates that mesh vertex counts changed between write and read, which should not happen.",
            n, dataNameStr, meshNameStr, buffer.size()));
      }
      for (std::size_t i = 0; i < n; ++i) {
        values[i] = buffer[i];
      }
    } else {
      // No buffer available - return zeros
      for (std::size_t i = 0; i < n; ++i) {
        values[i] = 0.0;
      }
    }
    break;
  }

  case impl::ParticipantImpl::DataMode::ScaledBuffer: {
    // Mode 3: Return buffered write data with scaling
    auto bufferIt = _impl->writeBuffers.find(key);
    if (bufferIt != _impl->writeBuffers.end() && !bufferIt->second.empty()) {
      const auto &buffer = bufferIt->second;
      // Strict validation: buffer size must match requested size
      if (buffer.size() != n) {
        throw precice::Error(precice::utils::format_or_error(
            "readData() was called with {} values for data '{}' on mesh '{}', "
            "but writeData() previously wrote {} values. "
            "The number of values in readData() must match writeData() for the same mesh/data pair. "
            "This typically indicates that mesh vertex counts changed between write and read, which should not happen.",
            n, dataNameStr, meshNameStr, buffer.size()));
      }

      if (!vectorMult.empty()) {
        // Element-wise multiplication with vector
        for (std::size_t i = 0; i < n; ++i) {
          values[i] = buffer[i] * vectorMult[i % vectorMult.size()];
        }
      } else {
        // Scalar multiplication
        for (std::size_t i = 0; i < n; ++i) {
          values[i] = buffer[i] * scalarMult;
        }
      }
    } else {
      // No buffer available - return zeros
      for (std::size_t i = 0; i < n; ++i) {
        values[i] = 0.0;
      }
    }
    break;
  }
  }
}

// Just-in-time mapping (experimental)

void Participant::writeAndMapData(
    ::precice::string_view /*meshName*/,
    ::precice::string_view /*dataName*/,
    precice::span<const double> coordinates,
    precice::span<const double> values)
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (!_impl->initialized) {
    throw precice::Error("initialize() has to be called before writeAndMapData().");
  }
  // coordinates are 3*N, values are N for scalar data
  if (coordinates.size() % 3 != 0) {
    throw precice::Error(precice::utils::format_or_error("writeAndMapData() was called with {} coordinates, but needs a multiple of 3 coordinates.", coordinates.size()));
  }
  if (values.size() != coordinates.size() / 3) {
    throw precice::Error(precice::utils::format_or_error("writeAndMapData() was called with {} coordinates and {} values, but needs {} values ({} / 3).", coordinates.size(), values.size(), coordinates.size() / 3, coordinates.size()));
  }
  // no-op
}

void Participant::mapAndReadData(
    ::precice::string_view /*meshName*/,
    ::precice::string_view /*dataName*/,
    precice::span<const double> coordinates,
    double /*relativeReadTime*/,
    precice::span<double> values) const
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (!_impl->initialized) {
    throw precice::Error("initialize() has to be called before mapAndReadData().");
  }
  if (coordinates.size() % 3 != 0) {
    throw precice::Error(precice::utils::format_or_error("mapAndReadData() was called with {} coordinates, but needs a multiple of 3 coordinates.", coordinates.size()));
  }
  if (values.size() != coordinates.size() / 3) {
    throw precice::Error(precice::utils::format_or_error("mapAndReadData() was called with {} coordinates and {} values, but needs {} values ({} / 3).", coordinates.size(), values.size(), coordinates.size() / 3, coordinates.size()));
  }
  // no-op
}

// Direct access

void Participant::setMeshAccessRegion(
    ::precice::string_view /*meshName*/,
    precice::span<const double> boundingBox) const
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (!_impl->initialized) {
    throw precice::Error("initialize() has to be called before setMeshAccessRegion().");
  }
  // bounding box: min/max for each dimension -> 3D: 6 entries
  if (boundingBox.size() != 6) {
    throw precice::Error(precice::utils::format_or_error("setMeshAccessRegion() was called with {} bounding box coordinates, but expects 6 coordinates (min and max for each dimension in 3D).", boundingBox.size()));
  }
  // no-op
}

void Participant::getMeshVertexIDsAndCoordinates(
    ::precice::string_view /*meshName*/,
    precice::span<VertexID> ids,
    precice::span<double>   coordinates) const
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (!_impl->initialized) {
    throw precice::Error("initialize() has to be called before getMeshVertexIDsAndCoordinates().");
  }
  // coordinates are 3 * ids.size()
  if (coordinates.size() != ids.size() * 3) {
    throw precice::Error(precice::utils::format_or_error("getMeshVertexIDsAndCoordinates() was called with {} vertices and {} coordinates, but needs {} coordinates ({} x 3).", ids.size(), coordinates.size(), ids.size() * 3, ids.size()));
  }
  // no-op
}

// Gradient data (experimental)

bool Participant::requiresGradientDataFor(::precice::string_view /*meshName*/,
                                          ::precice::string_view /*dataName*/) const
{
  return false;
}

void Participant::writeGradientData(
    ::precice::string_view /*meshName*/,
    ::precice::string_view /*dataName*/,
    precice::span<const VertexID> ids,
    precice::span<const double>   gradients)
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (!_impl->initialized) {
    throw precice::Error("initialize() has to be called before writeGradientData().");
  }
  // gradients expected as 3 components per id for vector gradient in 3D
  if (gradients.size() != ids.size() * 3) {
    throw precice::Error(precice::utils::format_or_error("writeGradientData() was called with {} vertices and {} gradient components, but needs {} components ({} x 3).", ids.size(), gradients.size(), ids.size() * 3, ids.size()));
  }
  // no-op
}

// Profiling

void Participant::startProfilingSection(::precice::string_view /*sectionName*/)
{
  // no-op
}

void Participant::stopLastProfilingSection()
{
  // no-op
}

// Missing Symbols

extern char const *const versionInformation __attribute__((visibility("default"))) = "precice-mock;revision=mock;features=MPI:OFF,PETSC:OFF,PYTHON:OFF"; // Tooling API

namespace tooling {

void printConfigReference(std::ostream &out, ConfigReferenceType reftype)
{
  out << "preCICE mock config reference; type=" << static_cast<int>(reftype) << "\n";
}

void checkConfiguration(const std::string &filename, const std::string &participant, int size)
{
  (void) filename;
  (void) participant;
  (void) size;
  std::cerr << "precice mock: checkConfiguration called\n";
}

} // namespace tooling

std::string getVersionInformation()
{
  return std::string("precice-mock;revision=mock;features=MPI:OFF,PETSC:OFF,PYTHON:OFF");
}

} // namespace precice
