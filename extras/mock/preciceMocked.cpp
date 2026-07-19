#include <algorithm>
#include <cmath>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <mutex>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <libxml/SAX2.h>

#include <precice/Exceptions.hpp>
#include <precice/Tooling.hpp>
#include <precice/Types.hpp>
#include <precice/impl/versions.hpp>
#include <precice/span.hpp>
#include "utils/fmt.hpp"

namespace {

std::string buildMockVersionInformation()
{
  // Lead with the version number to match real preCICE's format
  // (<version>;...). Tools such as the FMI runner parse the first character to
  // detect the major version, so a "precice-mock" prefix breaks them.
  std::string info = PRECICE_VERSION ";precice-mock";
#ifdef __VERSION__
  info += ";compiler=";
  info += __VERSION__;
#endif
#if defined(__cplusplus)
  info += ";cxx=";
  info += std::to_string(static_cast<long long>(__cplusplus));
#endif
  return info;
}

const std::string mockVersionInformation = buildMockVersionInformation();

// Tolerance used for time comparisons, matching preCICE's math::NUMERICAL_ZERO_DIFFERENCE.
// Like real preCICE (DoubleAggregator), the mock uses compensated summation for time
// accumulation, so the same tight tolerance is safe even with heavy subcycling.
constexpr double timeTolerance = 1e-14;

// Kahan (compensated) summation, mirroring preCICE's DoubleAggregator: keeps time
// accumulation errors at machine precision independent of the number of substeps.
struct KahanAccumulator {
  double sum  = 0.0;
  double comp = 0.0;

  void reset(double v = 0.0)
  {
    sum  = v;
    comp = 0.0;
  }
  void add(double v)
  {
    double y = v - comp;
    double t = sum + y;
    comp     = (t - sum) - y;
    sum      = t;
  }
};

std::string toStr(::precice::span<const char> sv)
{
  std::string s(sv.data(), sv.size());
  auto        pos = s.find('\0');
  if (pos != std::string::npos) {
    s.resize(pos);
  }
  return s;
}

} // namespace

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
  Performs basic checks and returns various basic data series for testing purposes.
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
    int         dimensions           = 3;
    bool        provided             = false;
    bool        received             = false;
    bool        requiresConnectivity = false;
    bool        apiAccess            = false;
  };

  struct DataInfo {
    std::string name;
    std::string meshName;
    int         dimensions = 1;
    bool        isRead     = false;
    bool        isWrite    = false;
  };

  struct ExchangeInfo {
    std::string dataName;
    std::string meshName;
    std::string fromParticipant;
    bool        initialize = false;
  };

  struct MappingInfo {
    std::string direction; // "read" or "write"
    std::string type;      // element name, e.g. "nearest-neighbor", "nearest-projection"
    std::string constraint;
    std::string fromMesh;
    std::string toMesh;
    bool        requiresGradient = false;
    bool        isJit            = false; // missing from/to attribute
  };

  struct ConfigData {
    std::map<std::string, MeshInfo> meshes;
    std::vector<DataInfo>           dataItems;
    std::vector<ExchangeInfo>       exchanges;
    std::vector<MappingInfo>        mappings;
    bool                            isImplicitCoupling             = false;
    bool                            hasConvergenceMeasure          = false;
    int                             minIterations                  = 1;
    int                             configMinIterations            = 1;
    int                             maxIterations                  = 1;
    int                             configMaxIterations            = 1;     // Stores the value from preCICE config
    int                             maxIterationsOverride          = 2;     // Default override (2), -1 use config, >0 override
    bool                            maxIterationsOverrideSpecified = false; // True if mock-config explicitly sets override
    double                          maxTime                        = -1.0;  // -1 means not set
    int                             maxTimeWindows                 = -1;    // -1 means not set
    double                          timeWindowSize                 = -1.0;  // -1 means not set
    // <time-window-size method="first-participant" />: the solver prescribes the
    // window size; every advance() completes a time window.
    bool timeWindowSizeFirstParticipant = false;
    bool allowsExperimental             = false;
    bool allowsRemeshing                = false;
  };

  // Internal runtime state separated from configuration
  struct RuntimeState {
    // Iteration/convergence
    int  currentIteration   = 0;
    bool iterationConverged = false;
    // Checkpoint flags
    bool writeCheckpointRequired = false;
    bool readCheckpointRequired  = false;
    // Participant lifecycle
    bool initialized         = false;
    bool finalized           = false;
    bool couplingOngoing     = false;
    bool timeWindowCompleted = false;
    // Timekeeping (compensated accumulation like preCICE's TimeHandler)
    double           timeWindowSize = 0.0; // effective runtime time-window size
    std::size_t      currentStep    = 0;
    KahanAccumulator windowStartTime;   // start time of the current window
    KahanAccumulator windowTime;        // time progressed inside the current window
    double           currentTime = 0.0; // windowStartTime + windowTime
    // Direct-access configuration
    std::map<std::string, std::vector<double>> meshAccessRegions;
    // Synthetic vertex coordinates handed out for direct-access meshes (whose
    // vertices come from the partner, not from local setMeshVertices()).
    std::map<std::string, std::vector<double>> directAccessCoords;
    // RNG seed
    uint32_t seed = 0;
    // Parsing/state flags
    bool configParsed     = false;
    bool mockConfigParsed = false;
    bool mockConfigExists = false;
    // Initial data
    bool initialDataRequired  = false;
    bool initialDataFulfilled = false;
    // Profiling
    std::size_t activeProfilingSections = 0;
    // Mesh lock: meshes are locked after initialize(), unlocked per-mesh by resetMesh()
    std::set<std::string> lockedMeshes;
    // Runtime containers
    // Write buffers keyed by "mesh:data". Reads look up the exact key first, then fall
    // back to the same data name on another mesh, then to the globally last-written data.
    std::map<std::string, std::vector<double>> writeBuffers;
    std::vector<double>                        lastWriteBuffer;
    std::map<std::string, std::size_t>         lastWriteSizes;
    std::map<std::string, std::size_t>         meshVertexCounts;
  };

  enum class DataMode {
    Random,      // Mode 1: Random data
    Buffer,      // Mode 2: Return written data buffer (default)
    ScaledBuffer // Mode 3: Return written data buffer scaled/multiplied
  };

  enum class LoggingMode {
    Mock,   // Mock-only: minimal logging, only essential mock behavior
    PrecICE // PrecICE-style: verbose logging mimicking real preCICE
  };

  struct MockDataConfig {
    std::string         dataName;
    std::string         meshName;
    DataMode            mode             = DataMode::Buffer; // Default to buffer mode
    double              scalarMultiplier = 1.0;
    double              randomLower      = 0.0;
    double              randomUpper      = 1.0;
    uint32_t            randomSeed       = 0; // 0 means use default seed from rank
    std::vector<double> vectorMultiplier;     // Empty means use scalar
  };

  struct MockConfig {
    std::map<std::string, MockDataConfig> dataConfigs; // Key: "meshName:dataName"
    DataMode                              defaultMode             = DataMode::Buffer;
    double                                defaultScalarMultiplier = 1.0;
    double                                defaultRandomLower      = 0.0;
    double                                defaultRandomUpper      = 1.0;
    uint32_t                              defaultRandomSeed       = 0;     // 0 means use default seed from rank
    std::vector<double>                   defaultVectorMultiplier;         // Empty means use scalar
    LoggingMode                           loggingMode = LoggingMode::Mock; // Default to mock mode
  };

  // Minimal internal state
  std::string                  name;
  std::string                  config;
  std::string                  mockConfigPath;
  int                          rank    = 0;
  bool                         primary = false;
  int                          size    = 1;
  void                        *comm    = nullptr;
  mutable std::recursive_mutex mtx;

  ConfigData   configData;
  RuntimeState runtimeState;

  MockConfig mockConfig;

  enum class DataType { Scalar = 1,
                        Vector = 3 };

  // SAX parsing state for config
  struct ConfigParseState {
    bool                            inParticipant = false;
    std::string                     participantName;
    std::map<std::string, DataType> dataType;                 // type marker: scalar or vector
    bool                            inCouplingScheme = false; // Track if we're inside a coupling-scheme element
  } configParseState;

  // Per-coupling-scheme accumulation. A config may contain several coupling schemes
  // (3+ participants); only schemes this participant is a member of are applied.
  struct SchemeParseState {
    bool                      isImplicit            = false;
    bool                      hasConvergenceMeasure = false;
    bool                      sawParticipants       = false; // scheme lists its participants
    bool                      isMember              = false; // this participant is listed
    int                       minIterations         = -1;
    int                       maxIterations         = -1;
    double                    maxTime               = -1.0;
    int                       maxTimeWindows        = -1;
    double                    timeWindowSize        = -1.0;
    bool                      twsFirstParticipant   = false;
    std::vector<ExchangeInfo> exchanges;
  } schemeParseState;

  // SAX parsing state for mock config
  struct MockParseState {
    bool                inMockedData = false;
    bool                inDefault    = false;
    std::string         currentMesh;
    std::string         currentData;
    DataMode            currentMode        = DataMode::Buffer;
    double              currentScalar      = 1.0;
    double              currentRandomLower = 0.0;
    double              currentRandomUpper = 1.0;
    uint32_t            currentRandomSeed  = 0;
    std::vector<double> currentVector;
  } mockParseState;

  // Configuration parsing methods
  void                     parseConfig();
  void                     parseMockConfig();
  void                     applyMaxIterationsOverride();
  void                     mergeCouplingScheme();
  std::string              logPrefix() const;
  std::string              formatCouplingState() const;
  std::vector<std::string> buildInitializeLogMessages() const;

  // Validation / query helpers
  bool isMeshUsed(const std::string &meshName) const;
  // Throws the real-preCICE "unknown mesh" error if this participant does not use the mesh.
  void requireMeshUsed(const std::string &meshName) const;
  // Time until the end of the current window (truncated by max-time), 0 once coupling
  // ended; mirrors BaseCouplingScheme::getNextTimeStepMaxSize.
  double nextTimeStepMaxSize() const;

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
    : name(toStr(participantName)),
      config(toStr(configurationFileName)),
      rank(solverProcessIndex),
      primary(rank == 0),
      size(solverProcessSize),
      comm(communicator)
{
  // Basic invariant checks for the mock participant
  if (name.empty()) {
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
  // Initialize seed for random data generation 0x9e3779b9u is the fractional part of the golden ratio scaled to 32 bits used as an arbitrary constant
  runtimeState.seed = static_cast<uint32_t>(rank) ^ 0x9e3779b9u;
}

std::string impl::ParticipantImpl::logPrefix() const
{
  return (mockConfig.loggingMode == LoggingMode::PrecICE) ? "---[precice]  " : "[precice-mock]  ";
}

std::string impl::ParticipantImpl::formatCouplingState() const
{
  std::ostringstream state;
  if (configData.isImplicitCoupling) {
    state << "it " << runtimeState.currentIteration;
    if (configData.minIterations > 0 && configData.maxIterations > 0) {
      state << " (min: " << configData.minIterations << ", max: " << configData.maxIterations << ")";
    } else if (configData.maxIterations > 0) {
      state << " (max: " << configData.maxIterations << ")";
    } else if (configData.minIterations > 0) {
      state << " (min: " << configData.minIterations << ")";
    }
    state << ", ";
  }

  const std::size_t timeWindowIndex = runtimeState.currentStep + 1;
  state << "time-window " << timeWindowIndex;
  if (configData.maxTimeWindows > 0) {
    state << " (max: " << configData.maxTimeWindows << ")";
  }
  state << ", t " << runtimeState.currentTime;
  if (configData.maxTime > 0) {
    state << " (max: " << configData.maxTime << ")";
  }
  if (configData.timeWindowSizeFirstParticipant) {
    state << ", Dt first-participant";
  } else {
    state << ", Dt " << runtimeState.timeWindowSize << ", max-dt " << (runtimeState.timeWindowSize - runtimeState.windowTime.sum);
  }
  return state.str();
}

std::vector<std::string> impl::ParticipantImpl::buildInitializeLogMessages() const
{
  std::vector<std::string> messages;
  const double             effectiveTimeWindowSize = (configData.timeWindowSize > 0) ? configData.timeWindowSize : 1.0;

  if (mockConfig.loggingMode == LoggingMode::PrecICE) {
    messages.push_back("This is preCICE version 3.3.0 (mock)");
    messages.push_back("Revision info: mock-implementation");
    messages.push_back("Build type: Release (mock)");
  }

  try {
    messages.push_back(std::string("Working directory \"") + std::filesystem::current_path().string() + "\"");
  } catch (const std::filesystem::filesystem_error &fse) {
    messages.push_back(std::string("Working directory unknown due to error \"") + fse.what() + "\"");
  }

  if (runtimeState.mockConfigExists) {
    messages.push_back(std::string("Configuring preCICE (mock) with configuration \"") +
                       std::filesystem::absolute(config).string() + "\" and mock configuration \"" +
                       std::filesystem::absolute(mockConfigPath).string() + "\"");
  } else {
    messages.push_back(std::string("Configuring preCICE (mock) with configuration \"") +
                       std::filesystem::absolute(config).string() +
                       "\" and default mock configuration. See the mock-README about how to adjust this.");
  }

  if (mockConfig.loggingMode == LoggingMode::PrecICE) {
    messages.push_back(std::string("I am participant \"") + name + "\"");
  } else {
    std::ostringstream summary;
    summary << "Initialized with participant \"" << name << "\"";
    messages.push_back(summary.str());

    std::ostringstream configSummary;
    configSummary << "Configuration: default-data-mode=";
    if (mockConfig.defaultMode == DataMode::Random) {
      configSummary << "random";
    } else if (mockConfig.defaultMode == DataMode::ScaledBuffer) {
      configSummary << "scaled";
    } else {
      configSummary << "buffer";
    }
    if (configData.isImplicitCoupling) {
      configSummary << " max-iterations=" << configData.maxIterations;
    }
    if (configData.timeWindowSizeFirstParticipant) {
      configSummary << " time-window-size=first-participant";
    } else {
      configSummary << " time-window-size=" << effectiveTimeWindowSize;
    }
    if (runtimeState.mockConfigExists) {
      configSummary << " mock-config=yes";
    }
    messages.push_back(configSummary.str());
  }

  if (mockConfig.loggingMode == LoggingMode::PrecICE) {
    messages.push_back("Setting up primary communication to coupling partner/s");
    messages.push_back("Primary ranks are connected");
    messages.push_back("Setting up preliminary secondary communication to coupling partner/s");

    for (const auto &meshEntry : configData.meshes) {
      const auto &meshInfo = meshEntry.second;
      if (!meshInfo.provided) {
        continue;
      }

      messages.push_back(std::string("Prepare partition for mesh ") + meshInfo.name);
      messages.push_back(std::string("Gather mesh ") + meshInfo.name);
      messages.push_back(std::string("Send global mesh ") + meshInfo.name);
    }

    messages.push_back(formatCouplingState());
  }

  return messages;
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

  // Handle precice-configuration root element. Its name contains a hyphen, not a
  // colon, so libxml2 reports it as a prefix-less element named "precice-configuration".
  if ((nsPrefix.empty() && elemName == "precice-configuration") ||
      (nsPrefix == "precice" && elemName == "configuration")) {
    impl->configData.allowsExperimental = attrIsTrue("experimental");
    impl->configData.allowsRemeshing    = attrIsTrue("allow-remeshing");
  }

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
  // Handle participant (the <participant> tags inside coupling-scheme:multi are
  // membership markers, handled in the coupling-scheme chain below)
  else if (elemName == "participant" && !impl->configParseState.inCouplingScheme) {
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
        if (elemName == "receive-mesh" || attrIsTrue("receive")) {
          it->second.received = true;
          // "direct-access" is the deprecated alias for "api-access" (still
          // accepted by real preCICE 3.x, removed in v4); honor both.
          if (attrIsTrue("api-access") || attrIsTrue("direct-access"))
            it->second.apiAccess = true;
        }
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
  // Handle mapping tags of ALL participants (namespace prefix "mapping"): the
  // connectivity requirement of a mesh may originate from the partner's mapping.
  if (nsPrefix == "mapping") {
    impl::ParticipantImpl::MappingInfo mi;
    mi.direction        = getAttr("direction");
    mi.type             = elemName;
    mi.constraint       = getAttr("constraint");
    mi.fromMesh         = getAttr("from");
    mi.toMesh           = getAttr("to");
    mi.requiresGradient = (elemName.find("gradient") != std::string::npos);
    mi.isJit            = (mi.direction == "write" && mi.fromMesh.empty()) ||
               (mi.direction == "read" && mi.toMesh.empty());

    // Derive connectivity ("full mesh") requirements, mirroring the real mappings'
    // MeshRequirement::FULL: projection-based mappings need connectivity on the input
    // mesh (consistent), the output mesh (conservative), or both (scaled-consistent);
    // scaled-consistent variants of all mappings need both.
    auto markFull = [&impl](const std::string &meshName) {
      if (meshName.empty())
        return;
      auto it = impl->configData.meshes.find(meshName);
      if (it != impl->configData.meshes.end())
        it->second.requiresConnectivity = true;
    };
    const bool scaled     = mi.constraint.find("scaled-consistent") != std::string::npos;
    const bool projective = (elemName == "nearest-projection" || elemName == "linear-cell-interpolation");
    if (projective) {
      if (mi.constraint == "conservative") {
        markFull(mi.toMesh);
      } else if (scaled) {
        markFull(mi.fromMesh);
        markFull(mi.toMesh);
      } else { // consistent (default)
        markFull(mi.fromMesh);
      }
    } else if (scaled) {
      markFull(mi.fromMesh);
      markFull(mi.toMesh);
    }

    if (impl->configParseState.inParticipant) {
      impl->configData.mappings.push_back(mi);
    }
  }
  // Handle coupling schemes - detect when we enter a coupling-scheme element.
  // Settings are accumulated per scheme and merged at the scheme's end element,
  // but only if this participant is a member of the scheme (see mergeCouplingScheme).
  // The "multi" scheme is always implicit and is used to couple more than two participants.
  if (elemName == "serial-implicit" || elemName == "parallel-implicit" ||
      elemName == "serial-explicit" || elemName == "parallel-explicit" ||
      elemName == "multi") {
    impl->configParseState.inCouplingScheme = true;
    impl->schemeParseState                  = SchemeParseState{};
    impl->schemeParseState.isImplicit       = (elemName == "serial-implicit" || elemName == "parallel-implicit" || elemName == "multi");
  } else if (impl->configParseState.inCouplingScheme && elemName == "participants") {
    impl->schemeParseState.sawParticipants = true;
    if (getAttr("first") == impl->name || getAttr("second") == impl->name) {
      impl->schemeParseState.isMember = true;
    }
  } else if (impl->configParseState.inCouplingScheme && elemName == "participant") {
    // membership marker of a coupling-scheme:multi
    impl->schemeParseState.sawParticipants = true;
    if (getAttr("name") == impl->name) {
      impl->schemeParseState.isMember = true;
    }
  } else if (impl->configParseState.inCouplingScheme && elemName == "max-iterations") {
    std::string maxIterStr = getAttr("value");
    if (!maxIterStr.empty()) {
      impl->schemeParseState.maxIterations = std::stoi(maxIterStr);
    }
  } else if (impl->configParseState.inCouplingScheme && elemName == "min-iterations") {
    std::string minIterStr = getAttr("value");
    if (!minIterStr.empty()) {
      impl->schemeParseState.minIterations = std::stoi(minIterStr);
    }
  } else if (impl->configParseState.inCouplingScheme &&
             (elemName == "relative-convergence-measure" || elemName == "absolute-convergence-measure" ||
              elemName == "residual-relative-convergence-measure" || elemName == "min-iteration-convergence-measure")) {
    impl->schemeParseState.hasConvergenceMeasure = true;
  } else if (impl->configParseState.inCouplingScheme && elemName == "exchange") {
    std::string dataName = getAttr("data");
    std::string meshName = getAttr("mesh");
    std::string fromPart = getAttr("from");
    bool        init     = attrIsTrue("initialize");
    if (!dataName.empty() && !fromPart.empty()) {
      impl::ParticipantImpl::ExchangeInfo ei;
      ei.dataName        = dataName;
      ei.meshName        = meshName;
      ei.fromParticipant = fromPart;
      ei.initialize      = init;
      impl->schemeParseState.exchanges.push_back(ei);
    }
  }
  // Parse max-time, max-time-windows, and time-window-size (only inside coupling-scheme)
  if (impl->configParseState.inCouplingScheme && elemName == "max-time") {
    std::string maxTimeStr = getAttr("value");
    if (!maxTimeStr.empty()) {
      impl->schemeParseState.maxTime = std::stod(maxTimeStr);
    }
  } else if (impl->configParseState.inCouplingScheme && elemName == "max-time-windows") {
    std::string maxTimeWindowsStr = getAttr("value");
    if (!maxTimeWindowsStr.empty()) {
      impl->schemeParseState.maxTimeWindows = std::stoi(maxTimeWindowsStr);
    }
  } else if (impl->configParseState.inCouplingScheme && elemName == "time-window-size") {
    std::string timeWindowSizeStr = getAttr("value");
    if (!timeWindowSizeStr.empty()) {
      impl->schemeParseState.timeWindowSize = std::stod(timeWindowSizeStr);
    }
    if (getAttr("method") == "first-participant") {
      impl->schemeParseState.twsFirstParticipant = true;
    }
  }
}

void impl::ParticipantImpl::onConfigEndElement(void *ctx, const xmlChar *localname, const xmlChar *prefix, const xmlChar *URI)
{
  auto       *impl = static_cast<impl::ParticipantImpl *>(ctx);
  std::string elemName(reinterpret_cast<const char *>(localname));
  if (elemName == "participant" && !impl->configParseState.inCouplingScheme) {
    impl->configParseState.inParticipant = false;
  } else if (elemName == "serial-implicit" || elemName == "parallel-implicit" ||
             elemName == "serial-explicit" || elemName == "parallel-explicit" ||
             elemName == "multi") {
    impl->mergeCouplingScheme();
    impl->configParseState.inCouplingScheme = false;
  }
}

void impl::ParticipantImpl::mergeCouplingScheme()
{
  auto &s = schemeParseState;
  // Schemes naming their participants only apply to members; schemes without a
  // participant list (simplified configs) apply to everyone.
  if (s.sawParticipants && !s.isMember) {
    s = SchemeParseState{};
    return;
  }
  if (s.isImplicit) {
    // If this participant is part of several schemes (compositional coupling), any
    // implicit membership requires checkpoint handling, like in real preCICE.
    configData.isImplicitCoupling = true;
    if (s.minIterations > 0) {
      configData.minIterations       = s.minIterations;
      configData.configMinIterations = s.minIterations;
    }
    if (s.maxIterations > 0) {
      configData.maxIterations       = s.maxIterations;
      configData.configMaxIterations = s.maxIterations;
    }
  }
  configData.hasConvergenceMeasure = configData.hasConvergenceMeasure || s.hasConvergenceMeasure;
  // A compositional simulation runs until every member scheme reached its end,
  // so keep the largest termination bounds.
  if (s.maxTime > 0) {
    configData.maxTime = std::max(configData.maxTime, s.maxTime);
  }
  if (s.maxTimeWindows > 0) {
    configData.maxTimeWindows = std::max(configData.maxTimeWindows, s.maxTimeWindows);
  }
  if (s.timeWindowSize > 0) {
    configData.timeWindowSize = s.timeWindowSize;
  }
  if (s.twsFirstParticipant) {
    configData.timeWindowSizeFirstParticipant = true;
  }
  configData.exchanges.insert(configData.exchanges.end(), s.exchanges.begin(), s.exchanges.end());
  s = SchemeParseState{};
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
  if (runtimeState.configParsed)
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
        // Writing data does not make a direct-access (received + api-access)
        // mesh a provided one: its vertices come from the partner via
        // getMeshVertexIDsAndCoordinates(), not from local setMeshVertices().
        if (dataInfo.isWrite && !(it->second.received && it->second.apiAccess))
          it->second.provided = true;
        if (dataInfo.isRead)
          it->second.received = true;
      }
    }

    // Set defaults for implicit coupling
    if (configData.isImplicitCoupling) {
      // Only print warning if no convergence measure was found
      if (!configData.hasConvergenceMeasure) {
        std::cout << precice::utils::format_or_error(
                         "{}WARNING: No convergence measures were defined for an implicit coupling scheme. "
                         "It will always iterate the maximum amount iterations, which is {}. "
                         "You may want to add a convergence measure in your <coupling-scheme:.../> in your configuration.",
                         logPrefix(), configData.maxIterations)
                  << std::endl;
      }
    }

    // Validate termination configuration
    if (configData.maxTime < 0 && configData.maxTimeWindows < 0) {
      throw precice::Error(precice::utils::format_or_error(
          "No termination configuration found in preCICE config '{}'. "
          "Please add <max-time value=\"...\"/> or <max-time-windows value=\"...\"/> "
          "to your <coupling-scheme:.../> to prevent infinite execution.",
          config));
    }

    applyMaxIterationsOverride();

    if (configData.maxIterations > 0 && configData.minIterations > configData.maxIterations) {
      configData.minIterations = configData.maxIterations;
    }

    // Determine if initial data is required (exchange with initialize="true" from this participant)
    for (const auto &ex : configData.exchanges) {
      if (ex.initialize && ex.fromParticipant == name) {
        runtimeState.initialDataRequired = true;
        break;
      }
    }

    runtimeState.configParsed = true;
  } catch (const precice::Error &) {
    throw;
  } catch (const std::exception &e) {
    throw precice::Error(precice::utils::format_or_error(
        "{}Warning: Could not parse configuration file '{}': {}. Using limited mock behavior.",
        logPrefix(), config, e.what()));
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

  // Handle logging-mode at root level
  if (elemName == "logging-mode") {
    std::string modeStr = attrs["mode"];
    if (modeStr == "precice") {
      impl->mockConfig.loggingMode = LoggingMode::PrecICE;
    } else if (modeStr == "mock") {
      impl->mockConfig.loggingMode = LoggingMode::Mock;
    }
    return;
  }
  // Handle max-iterations-override at root level
  if (elemName == "max-iterations-override") {
    std::string valueStr = attrs["value"];
    if (!valueStr.empty()) {
      int override = std::stoi(valueStr);
      if (override == -1 || override > 0) {
        impl->configData.maxIterationsOverride          = override;
        impl->configData.maxIterationsOverrideSpecified = true;
      }
    }
  } else if (elemName == "mocked-data-default") {
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
    impl->mockParseState.currentRandomLower = impl->mockConfig.defaultRandomLower;
    impl->mockParseState.currentRandomUpper = impl->mockConfig.defaultRandomUpper;
    impl->mockParseState.currentRandomSeed  = impl->mockConfig.defaultRandomSeed;

    std::string lowerStr = attrs["lower"];
    std::string upperStr = attrs["upper"];
    if (!lowerStr.empty()) {
      impl->mockParseState.currentRandomLower = std::stod(lowerStr);
    }
    if (!upperStr.empty()) {
      impl->mockParseState.currentRandomUpper = std::stod(upperStr);
    }
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
    impl->mockParseState.currentRandomLower = impl->mockConfig.defaultRandomLower;
    impl->mockParseState.currentRandomUpper = impl->mockConfig.defaultRandomUpper;
    impl->mockParseState.currentRandomSeed  = impl->mockConfig.defaultRandomSeed;

    std::string lowerStr = attrs["lower"];
    std::string upperStr = attrs["upper"];
    if (!lowerStr.empty()) {
      impl->mockParseState.currentRandomLower = std::stod(lowerStr);
    }
    if (!upperStr.empty()) {
      impl->mockParseState.currentRandomUpper = std::stod(upperStr);
    }
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
    } else if (elemName == "bounds") {
      std::string lowerStr = attrs["lower"];
      std::string upperStr = attrs["upper"];
      if (!lowerStr.empty()) {
        impl->mockParseState.currentRandomLower = std::stod(lowerStr);
      }
      if (!upperStr.empty()) {
        impl->mockParseState.currentRandomUpper = std::stod(upperStr);
      }
    } else if (elemName == "seed") {
      std::string seedStr = attrs["value"];
      if (!seedStr.empty()) {
        impl->mockParseState.currentRandomSeed = static_cast<uint32_t>(std::stoul(seedStr));
      }
    }
  }
}

void impl::ParticipantImpl::onMockEndElement(void *ctx, const xmlChar *localname, const xmlChar *prefix, const xmlChar *URI)
{
  auto       *impl = static_cast<impl::ParticipantImpl *>(ctx);
  std::string elemName(reinterpret_cast<const char *>(localname));
  if (elemName == "mocked-data-default" && impl->mockParseState.inDefault) {
    // Apply default multipliers
    impl->mockConfig.defaultScalarMultiplier = impl->mockParseState.currentScalar;
    impl->mockConfig.defaultVectorMultiplier = impl->mockParseState.currentVector;
    impl->mockConfig.defaultRandomLower      = impl->mockParseState.currentRandomLower;
    impl->mockConfig.defaultRandomUpper      = impl->mockParseState.currentRandomUpper;
    impl->mockConfig.defaultRandomSeed       = impl->mockParseState.currentRandomSeed;
    if (impl->mockConfig.defaultMode == DataMode::Random) {
      if (!std::isfinite(impl->mockConfig.defaultRandomLower) || !std::isfinite(impl->mockConfig.defaultRandomUpper)) {
        throw precice::Error("mock-config default random bounds must be finite.");
      }
      if (impl->mockConfig.defaultRandomUpper <= impl->mockConfig.defaultRandomLower) {
        throw precice::Error(precice::utils::format_or_error(
            "mock-config default random bounds are invalid: upper={} must be greater than lower={}.",
            impl->mockConfig.defaultRandomUpper, impl->mockConfig.defaultRandomLower));
      }
    }
    if (impl->mockConfig.defaultMode == DataMode::ScaledBuffer) {
      bool allOnes = !impl->mockParseState.currentVector.empty();
      if (allOnes) {
        for (double v : impl->mockParseState.currentVector) {
          if (v != 1.0) {
            allOnes = false;
            break;
          }
        }
      }
      if ((impl->mockParseState.currentVector.empty() && impl->mockParseState.currentScalar == 1.0) || allOnes) {
        std::cout << impl->logPrefix() << "WARNING: mock-config default uses mode 'scaled' but the multiplier is effectively 1 (no scaling)." << std::endl;
      }
    }
    if (impl->mockConfig.defaultMode == DataMode::Buffer &&
        (impl->mockParseState.currentScalar != 1.0 || !impl->mockParseState.currentVector.empty())) {
      std::cout << impl->logPrefix() << "WARNING: mock-config default uses mode 'buffer' but also specifies scalar/vector multipliers; multipliers are ignored for buffer mode." << std::endl;
    }
    impl->mockParseState.inDefault = false;
  } else if (elemName == "mocked-data" && impl->mockParseState.inMockedData) {
    if (!impl->mockParseState.currentMesh.empty() && !impl->mockParseState.currentData.empty()) {
      std::string    key = impl->mockParseState.currentMesh + ":" + impl->mockParseState.currentData;
      MockDataConfig config;
      config.meshName         = impl->mockParseState.currentMesh;
      config.dataName         = impl->mockParseState.currentData;
      config.mode             = impl->mockParseState.currentMode;
      config.scalarMultiplier = impl->mockParseState.currentScalar;
      config.randomLower      = impl->mockParseState.currentRandomLower;
      config.randomUpper      = impl->mockParseState.currentRandomUpper;
      config.randomSeed       = impl->mockParseState.currentRandomSeed;
      config.vectorMultiplier = impl->mockParseState.currentVector;
      if (config.mode == DataMode::Random) {
        if (!std::isfinite(config.randomLower) || !std::isfinite(config.randomUpper)) {
          throw precice::Error(precice::utils::format_or_error(
              "mock-config for mesh '{}' data '{}' uses non-finite random bounds (lower={}, upper={}).",
              config.meshName, config.dataName, config.randomLower, config.randomUpper));
        }
        if (config.randomUpper <= config.randomLower) {
          throw precice::Error(precice::utils::format_or_error(
              "mock-config for mesh '{}' data '{}' uses invalid random bounds: upper={} must be greater than lower={}.",
              config.meshName, config.dataName, config.randomUpper, config.randomLower));
        }
      }
      if (config.mode == DataMode::ScaledBuffer) {
        bool allOnes = !config.vectorMultiplier.empty();
        if (allOnes) {
          for (double v : config.vectorMultiplier) {
            if (v != 1.0) {
              allOnes = false;
              break;
            }
          }
        }
        if ((config.vectorMultiplier.empty() && config.scalarMultiplier == 1.0) || allOnes) {
          std::cout << impl->logPrefix() << "WARNING: mock-config for mesh '" << config.meshName << "' data '" << config.dataName
                    << "' uses mode 'scaled' but the multiplier is effectively 1 (no scaling)." << std::endl;
        }
      }
      if (config.mode == DataMode::Buffer &&
          (config.scalarMultiplier != 1.0 || !config.vectorMultiplier.empty())) {
        std::cout << impl->logPrefix() << "WARNING: mock-config for mesh '" << config.meshName << "' data '" << config.dataName
                  << "' uses mode 'buffer' but also specifies scalar/vector multipliers; multipliers are ignored for buffer mode." << std::endl;
      }
      impl->mockConfig.dataConfigs[key] = config;
    }
    impl->mockParseState.inMockedData = false;
  }
}

void impl::ParticipantImpl::parseMockConfig()
{
  if (runtimeState.mockConfigParsed)
    return;

  try {
    // Derive mock-config path from precice config path
    mockConfigPath   = config;
    size_t lastSlash = mockConfigPath.find_last_of("/\\");
    size_t lastDot   = mockConfigPath.find_last_of('.');

    if (lastSlash != std::string::npos) {
      if (lastDot != std::string::npos && lastDot > lastSlash) {
        mockConfigPath = mockConfigPath.substr(0, lastSlash + 1) + "precice-mock-config.xml";
      } else {
        mockConfigPath = mockConfigPath + "/precice-mock-config.xml";
      }
    } else {
      // No directory, same directory as executable
      mockConfigPath = "mock-config.xml";
    }

    // Check if mock config file exists
    std::ifstream mockFile(mockConfigPath);
    if (!mockFile.good()) {
      // No mock config - use defaults
      mockConfig.defaultMode        = DataMode::Buffer;
      runtimeState.mockConfigParsed = true;
      runtimeState.mockConfigExists = false;
      return;
    }
    mockFile.close();
    runtimeState.mockConfigExists = true;

    // Parse mock configuration file
    parseXMLFile(mockConfigPath, onMockStartElement, onMockEndElement);

    applyMaxIterationsOverride();
    runtimeState.mockConfigParsed = true;

  } catch (const std::exception &e) {
    throw precice::Error(precice::utils::format_or_error(
        "{}Warning: Could not parse mock configuration file: {}. Using default buffer mode.",
        logPrefix(), e.what()));
  }
}

void impl::ParticipantImpl::applyMaxIterationsOverride()
{
  int override = configData.maxIterationsOverride;

  // If mock-config did not explicitly specify an override and the preCICE
  // config sets max-iterations to 1, keep the config value (do not override).
  if (!configData.maxIterationsOverrideSpecified && configData.configMaxIterations == 1) {
    configData.maxIterations = configData.configMaxIterations;
    return;
  }

  if (override == 0) {
    // No override specified; keep config value
    return;
  }

  if (override == -1) {
    // Respect preCICE config value
    configData.maxIterations = configData.configMaxIterations;
    return;
  }

  if (override > 0) {
    if (configData.configMaxIterations > 0 && override > configData.configMaxIterations) {
      std::cerr << logPrefix() << "WARNING: max-iterations-override (" << override
                << ") is higher than max-iterations from preCICE config ("
                << configData.configMaxIterations << "). Using override value." << std::endl;
    }
    configData.maxIterations = override;
  }
}

bool impl::ParticipantImpl::isMeshUsed(const std::string &meshName) const
{
  auto it = configData.meshes.find(meshName);
  return it != configData.meshes.end() && (it->second.provided || it->second.received);
}

void impl::ParticipantImpl::requireMeshUsed(const std::string &meshName) const
{
  if (isMeshUsed(meshName)) {
    return;
  }
  // Mirrors the real library's PRECICE_VALIDATE_MESH_NAME error including the hint.
  std::vector<std::string> used;
  for (const auto &m : configData.meshes) {
    if (m.second.provided || m.second.received) {
      used.push_back(m.first);
    }
  }
  std::string hint;
  if (used.size() == 1) {
    hint = " This participant only knows mesh \"" + used.front() + "\".";
  } else if (!used.empty()) {
    hint = " Available meshes are: ";
    for (std::size_t i = 0; i < used.size(); ++i) {
      if (i > 0)
        hint += ", ";
      hint += used[i];
    }
  }
  throw precice::Error(precice::utils::format_or_error(
      "The mesh named \"{}\" is unknown to preCICE.{}", meshName, hint));
}

double impl::ParticipantImpl::nextTimeStepMaxSize() const
{
  if (!runtimeState.couplingOngoing) {
    return 0.0;
  }
  double maxDt;
  if (configData.timeWindowSizeFirstParticipant) {
    // The solver prescribes the window size: only max-time limits the step.
    maxDt = (configData.maxTime > 0) ? configData.maxTime - runtimeState.currentTime
                                     : std::numeric_limits<double>::max();
  } else {
    maxDt = runtimeState.timeWindowSize - runtimeState.windowTime.sum;
    if (configData.maxTime > 0) {
      maxDt = std::min(maxDt, configData.maxTime - runtimeState.currentTime);
    }
  }
  return maxDt;
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
  _impl->parseMockConfig();
  _impl->parseConfig();
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
  if (communicator == nullptr) {
    throw precice::Error("Passing \"nullptr\" as \"communicator\" to Participant constructor is not allowed. "
                         "Please use the Participant constructor without the \"communicator\" argument, if you don't want to pass an MPI communicator.");
  }
  _impl->parseMockConfig();
  _impl->parseConfig();
}

Participant::~Participant()
{
  if (!_impl->runtimeState.finalized) {
    finalize();
  }
}

// Steering methods

void Participant::initialize()
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (_impl->runtimeState.finalized) {
    throw precice::Error("initialize() cannot be called after finalize().");
  }
  if (_impl->runtimeState.initialized) {
    throw precice::Error("initialize() may only be called once.");
  }
  if (_impl->runtimeState.activeProfilingSections > 0) {
    throw precice::Error("There are unstopped user defined events. Please stop them using stopLastProfilingSection() before calling initialize().");
  }

  if (_impl->runtimeState.initialDataRequired && !_impl->runtimeState.initialDataFulfilled) {
    throw precice::Error("Initial data has to be written to preCICE before calling initialize(). "
                         "After defining your mesh, call requiresInitialData() to check if the participant is required to write initial data using the writeData() function.");
  }

  // Mesh preparation checks for provided meshes
  for (const auto &meshEntry : _impl->configData.meshes) {
    const auto &meshInfo = meshEntry.second;
    if (!meshInfo.provided) {
      continue;
    }

    std::size_t count = 0;
    auto        it    = _impl->runtimeState.meshVertexCounts.find(meshInfo.name);
    if (it != _impl->runtimeState.meshVertexCounts.end()) {
      count = it->second;
    }
    if (count == 0) {
      throw precice::Error(precice::utils::format_or_error(
          "The provided mesh \"{}\" is empty. Please set the mesh using setMeshVertex()/setMeshVertices() prior to calling initialize().",
          meshInfo.name));
    }
  }

  _impl->runtimeState.initialized         = true;
  _impl->runtimeState.couplingOngoing     = true;
  _impl->runtimeState.timeWindowCompleted = false;
  // With method="first-participant" the solver prescribes the window size; the stored
  // value is unused then (see nextTimeStepMaxSize()).
  _impl->runtimeState.timeWindowSize = (_impl->configData.timeWindowSize > 0) ? _impl->configData.timeWindowSize : 1.0;
  _impl->runtimeState.currentStep    = 0;
  _impl->runtimeState.windowStartTime.reset();
  _impl->runtimeState.windowTime.reset();
  _impl->runtimeState.currentTime = 0.0;

  // Initialize implicit coupling state
  if (_impl->configData.isImplicitCoupling) {
    _impl->runtimeState.currentIteration        = 1;
    _impl->runtimeState.iterationConverged      = false;
    _impl->runtimeState.writeCheckpointRequired = true;
    _impl->runtimeState.readCheckpointRequired  = false;
  }

  // Validate that termination criteria are provided
  if (_impl->configData.maxTime <= 0 && _impl->configData.maxTimeWindows <= 0) {
    throw precice::Error(
        "The preCICE configuration must specify either max-time or max-time-windows in the coupling-scheme \n"
        "to prevent infinite execution. Please add one of these attributes:\n"
        "  <max-time value=\"...\"/>\n"
        "  <max-time-windows value=\"...\"/>");
  }

  // Lock all meshes after initialization
  for (const auto &meshEntry : _impl->configData.meshes) {
    _impl->runtimeState.lockedMeshes.insert(meshEntry.first);
  }

  // Print startup banner (only on rank 0) at the end of the setup block.
  if (_impl->primary) {
    const std::string prefix = _impl->logPrefix();
    for (const auto &message : _impl->buildInitializeLogMessages()) {
      std::cout << prefix << message << std::endl;
    }
  }
}

void Participant::advance(double computedTimeStepSize)
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (_impl->runtimeState.activeProfilingSections > 0) {
    throw precice::Error("There are unstopped user defined events. Please stop them using stopLastProfilingSection() before calling advance().");
  }
  if (_impl->runtimeState.finalized) {
    throw precice::Error("advance() cannot be called after finalize().");
  }
  if (!_impl->runtimeState.initialized) {
    throw precice::Error("initialize() has to be called before advance().");
  }
  if (!isCouplingOngoing()) {
    throw precice::Error("advance() cannot be called when isCouplingOngoing() returns false.");
  }
  if (!std::isfinite(computedTimeStepSize)) {
    throw precice::Error("advance() cannot be called with an infinite time step size.");
  }
  if (std::abs(computedTimeStepSize) < 1e-15) {
    throw precice::Error("advance() cannot be called with a time step size of 0.");
  }
  if (computedTimeStepSize < 0.0) {
    throw precice::Error(precice::utils::format_or_error("advance() cannot be called with a negative time step size {}.", computedTimeStepSize));
  }
  if (_impl->configData.timeWindowSizeFirstParticipant &&
      computedTimeStepSize >= std::numeric_limits<double>::max()) {
    throw precice::Error(
        "advance() was called with the max value of double which is not allowed. "
        "As this participant prescribes the time-window size using <time-window-size method=\"first-participant\" />, directly using getMaxTimeStepSize() is not permitted. "
        "Make sure to pass your own desired time-step size or use the recommended limiting \"dt = min(solver_dt, getMaxTimeStepSize())\".");
  }

  const double remaining = _impl->nextTimeStepMaxSize();
  // Uses a tolerance like preCICE's math::greaterEquals to permit dt == remaining.
  if (computedTimeStepSize - remaining > timeTolerance) {
    throw precice::Error(precice::utils::format_or_error(
        "The time step size given to preCICE in \"advance\" {} exceeds the maximum allowed time step size {} "
        "in the remaining of this time window. "
        "Did you restrict your time step size, \"dt = min(preciceDt, solverDt)\"? "
        "For more information, consult the adapter example in the preCICE documentation.",
        computedTimeStepSize, remaining));
  }

  // Enforce that required checkpoint actions were acknowledged via the requires*()
  // API before advancing, mirroring BaseCouplingScheme::checkCompletenessRequiredActions.
  if (_impl->configData.isImplicitCoupling) {
    if (_impl->runtimeState.readCheckpointRequired) {
      throw precice::Error("The iteration checkpoint wasn't read. Did you forget to call \"requiresReadingCheckpoint()\"?");
    }
    if (_impl->runtimeState.writeCheckpointRequired) {
      throw precice::Error("The iteration checkpoint wasn't written. Did you forget to call \"requiresWritingCheckpoint()\"?");
    }
  }

  // Re-lock all meshes (unlocked by resetMesh during this window)
  for (const auto &meshEntry : _impl->configData.meshes) {
    _impl->runtimeState.lockedMeshes.insert(meshEntry.first);
  }

  _impl->runtimeState.timeWindowCompleted = false;
  // isTimeWindowComplete() reports whether the *last* advance completed a window,
  // so the implicit-coupling flag has to be cleared here as well (it would otherwise
  // stay true throughout subcycled steps of the following window).
  _impl->runtimeState.iterationConverged = false;
  auto printCouplingState                = [&]() {
    if (_impl->primary && _impl->mockConfig.loggingMode == impl::ParticipantImpl::LoggingMode::PrecICE) {
      std::cout << _impl->logPrefix() << _impl->formatCouplingState() << std::endl;
    }
  };

  auto handleTimeWindowCompletion = [&]() {
    _impl->runtimeState.timeWindowCompleted = true;
    // Advance the window start by the actual progress. For regular windows this equals
    // timeWindowSize; for a truncated final window (when max-time is reached and is not a
    // multiple of the window size) this is smaller, matching TimeHandler::completeTimeWindow.
    _impl->runtimeState.windowStartTime.add(_impl->runtimeState.windowTime.sum);
    _impl->runtimeState.windowTime.reset();
    _impl->runtimeState.currentTime = _impl->runtimeState.windowStartTime.sum;
    _impl->runtimeState.currentStep += 1;

    bool shouldTerminate = false;
    if (_impl->configData.maxTimeWindows > 0 && _impl->runtimeState.currentStep >= static_cast<std::size_t>(_impl->configData.maxTimeWindows)) {
      shouldTerminate = true;
    }
    // Reached max-time: compare with tolerance like TimeHandler::reachedEnd.
    if (_impl->configData.maxTime > 0 && _impl->configData.maxTime - _impl->runtimeState.currentTime <= timeTolerance) {
      shouldTerminate = true;
    }

    if (shouldTerminate) {
      _impl->runtimeState.couplingOngoing = false;
    }
  };

  // A time window ends when no time remains until its end (within tolerance), mirroring
  // TimeHandler::reachedEndOfWindow which also accounts for truncation by max-time.
  // When the solver prescribes the window size (first-participant), every step ends one.
  auto reachedEndOfWindow = [&]() {
    if (_impl->configData.timeWindowSizeFirstParticipant) {
      return true;
    }
    return _impl->nextTimeStepMaxSize() <= timeTolerance;
  };

  // Handle implicit coupling iterations
  if (_impl->configData.isImplicitCoupling) {

    _impl->runtimeState.windowTime.add(computedTimeStepSize);
    _impl->runtimeState.currentTime = _impl->runtimeState.windowStartTime.sum + _impl->runtimeState.windowTime.sum;

    if (reachedEndOfWindow()) {
      if (_impl->runtimeState.currentIteration < _impl->configData.maxIterations) {
        // Start the next implicit coupling iteration within the same time window.
        _impl->runtimeState.currentIteration++;
        _impl->runtimeState.iterationConverged     = false;
        _impl->runtimeState.readCheckpointRequired = true;
        _impl->runtimeState.windowTime.reset();
        _impl->runtimeState.currentTime = _impl->runtimeState.windowStartTime.sum;
      } else {
        // The last implicit iteration completed the time window.
        _impl->runtimeState.iterationConverged      = true;
        _impl->runtimeState.currentIteration        = 1;
        _impl->runtimeState.writeCheckpointRequired = true;
        handleTimeWindowCompletion();
      }
    }
  } else {
    // Explicit coupling: always advance
    _impl->runtimeState.windowTime.add(computedTimeStepSize);
    _impl->runtimeState.currentTime = _impl->runtimeState.windowStartTime.sum + _impl->runtimeState.windowTime.sum;

    if (reachedEndOfWindow()) {
      handleTimeWindowCompletion();
    }
  }

  if (_impl->primary && _impl->mockConfig.loggingMode == impl::ParticipantImpl::LoggingMode::PrecICE) {
    const std::string prefix = _impl->logPrefix();
    if (_impl->runtimeState.timeWindowCompleted) {
      std::cout << prefix << "Time window completed" << std::endl;
    }

    if (!_impl->runtimeState.couplingOngoing) {
      std::cout << prefix << "Reached end at: final time-window: " << _impl->runtimeState.currentStep << ", final time: "
                << _impl->runtimeState.currentTime << std::endl;
    } else {
      printCouplingState();
    }
  }
}

void Participant::finalize()
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (_impl->runtimeState.finalized) {
    throw precice::Error("finalize() may only be called once.");
  }

  // Print finalize message (only on rank 0)
  if (_impl->primary) {
    const std::string prefix = _impl->logPrefix();
    if (_impl->mockConfig.loggingMode == impl::ParticipantImpl::LoggingMode::PrecICE) {
      std::cout << prefix << "Finalizing preCICE (mock)" << std::endl;
    } else {
      std::cout << prefix << "Finalized" << std::endl;
    }
  }

  _impl->runtimeState.initialized     = false;
  _impl->runtimeState.couplingOngoing = false;
  _impl->runtimeState.finalized       = true;
}

// Implicit coupling

bool Participant::requiresWritingCheckpoint()
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (_impl->runtimeState.finalized) {
    throw precice::Error("requiresWritingCheckpoint() cannot be called after finalize().");
  }
  if (!_impl->runtimeState.initialized) {
    throw precice::Error("initialize() has to be called before requiresWritingCheckpoint().");
  }

  // For implicit coupling, require checkpoint at start of iteration
  if (_impl->configData.isImplicitCoupling) {
    const bool required                         = _impl->runtimeState.writeCheckpointRequired;
    _impl->runtimeState.writeCheckpointRequired = false;
    return required;
  }

  return false;
}

bool Participant::requiresReadingCheckpoint()
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (_impl->runtimeState.finalized) {
    throw precice::Error("requiresReadingCheckpoint() cannot be called after finalize().");
  }
  if (!_impl->runtimeState.initialized) {
    throw precice::Error("initialize() has to be called before requiresReadingCheckpoint().");
  }

  // For implicit coupling, require reading checkpoint when not converged
  if (_impl->configData.isImplicitCoupling) {
    const bool required                        = _impl->runtimeState.readCheckpointRequired;
    _impl->runtimeState.readCheckpointRequired = false;
    return required;
  }

  return false;
}

// Status queries

int Participant::getMeshDimensions(::precice::string_view meshName) const
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);

  std::string meshNameStr = toStr(meshName);

  // Check if config was parsed and mesh exists
  if (_impl->runtimeState.configParsed) {
    _impl->requireMeshUsed(meshNameStr);
    return _impl->configData.meshes.at(meshNameStr).dimensions;
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

  std::string meshNameStr = toStr(meshName);
  std::string dataNameStr = toStr(dataName);

  // Check if config was parsed and data exists
  if (_impl->runtimeState.configParsed) {
    _impl->requireMeshUsed(meshNameStr);
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
  if (_impl->runtimeState.finalized) {
    throw precice::Error("isCouplingOngoing() cannot be called after finalize().");
  }
  if (!_impl->runtimeState.initialized) {
    throw precice::Error("initialize() has to be called before isCouplingOngoing() can be evaluated.");
  }
  return _impl->runtimeState.couplingOngoing;
}

bool Participant::isTimeWindowComplete() const
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (!_impl->runtimeState.initialized) {
    throw precice::Error("initialize() has to be called before isTimeWindowComplete().");
  }
  if (_impl->runtimeState.finalized) {
    throw precice::Error("isTimeWindowComplete() cannot be called after finalize().");
  }

  if (_impl->configData.isImplicitCoupling) {
    return _impl->runtimeState.iterationConverged;
  }

  return _impl->runtimeState.timeWindowCompleted;
}

double Participant::getMaxTimeStepSize() const
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (_impl->runtimeState.finalized) {
    throw precice::Error("getMaxTimeStepSize() cannot be called after finalize().");
  }
  if (!_impl->runtimeState.initialized) {
    throw precice::Error("initialize() has to be called before getMaxTimeStepSize() can be evaluated.");
  }
  // Matches BaseCouplingScheme::getNextTimeStepMaxSize: zero once coupling has ended, otherwise
  // the time until the end of the window, truncated by max-time. With
  // <time-window-size method="first-participant" /> the solver prescribes the step,
  // so only max-time limits it.
  return _impl->nextTimeStepMaxSize();
}

// Mesh access

bool Participant::requiresMeshConnectivityFor(::precice::string_view meshName) const
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  std::string                           meshNameStr = toStr(meshName);
  _impl->requireMeshUsed(meshNameStr);
  return _impl->configData.meshes.at(meshNameStr).requiresConnectivity;
}

void Participant::resetMesh(::precice::string_view meshName)
{
  std::lock_guard<std::recursive_mutex> lkm(_impl->mtx);
  if (!_impl->configData.allowsExperimental) {
    throw precice::Error("You called the API function \"resetMesh\", which is part of the experimental API. "
                         "You may unlock the full API by specifying <precice-configuration experimental=\"true\" ... > in the configuration. Please be aware that experimental features may change in any future version (even minor or bugfix).");
  }
  if (!_impl->configData.allowsRemeshing) {
    throw precice::Error("Cannot reset meshes. This feature needs to be enabled using <precice-configuration experimental=\"1\" allow-remeshing=\"1\">.");
  }
  if (_impl->runtimeState.finalized) {
    throw precice::Error("resetMesh() cannot be called after finalize().");
  }
  if (!_impl->runtimeState.initialized) {
    throw precice::Error("initialize() has to be called before resetMesh().");
  }
  if (!_impl->runtimeState.couplingOngoing) {
    throw precice::Error("Cannot remesh after the last time window has been completed.");
  }
  if (_impl->configData.isImplicitCoupling && !_impl->runtimeState.iterationConverged) {
    throw precice::Error("Cannot remesh while subcycling or iterating. Remeshing is only allowed when the time window is completed.");
  }
  std::string meshNameStr = toStr(meshName);
  _impl->requireMeshUsed(meshNameStr);
  _impl->runtimeState.lockedMeshes.erase(meshNameStr);
  _impl->runtimeState.meshVertexCounts[meshNameStr] = 0;
}

VertexID Participant::setMeshVertex(
    ::precice::string_view      meshName,
    precice::span<const double> position)
{
  std::lock_guard<std::recursive_mutex> lkm(_impl->mtx);
  std::string                           meshNameStr = toStr(meshName);
  _impl->requireMeshUsed(meshNameStr);
  if (_impl->runtimeState.lockedMeshes.count(meshNameStr)) {
    throw precice::Error(precice::utils::format_or_error(
        "This participant attempted to modify the Mesh \"{}\" while locked. "
        "Mesh modification is only allowed before calling initialize().",
        meshNameStr));
  }
  const int expectedDims = _impl->configData.meshes.at(meshNameStr).dimensions;
  if (position.size() != static_cast<std::size_t>(expectedDims)) {
    throw precice::Error(precice::utils::format_or_error("setMeshVertex() was called with {} coordinates, but expects {} coordinates for this mesh.", position.size(), expectedDims));
  }
  auto &count = _impl->runtimeState.meshVertexCounts[meshNameStr];
  ++count;
  return VertexID(count - 1);
}

int Participant::getMeshVertexSize(::precice::string_view meshName) const
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  std::string                           meshNameStr = toStr(meshName);
  _impl->requireMeshUsed(meshNameStr);
  // Received meshes only carry data after initialize(), like in the real library.
  const auto &meshInfo = _impl->configData.meshes.at(meshNameStr);
  if (!_impl->runtimeState.initialized && meshInfo.received && !meshInfo.provided) {
    throw precice::Error(precice::utils::format_or_error(
        "initialize() has to be called before accessing data of the received mesh \"{}\" on participant \"{}\".",
        meshNameStr, _impl->name));
  }
  auto it = _impl->runtimeState.meshVertexCounts.find(meshNameStr);
  if (it != _impl->runtimeState.meshVertexCounts.end()) {
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
  std::string                           meshNameStr = toStr(meshName);
  _impl->requireMeshUsed(meshNameStr);
  if (_impl->runtimeState.lockedMeshes.count(meshNameStr)) {
    throw precice::Error(precice::utils::format_or_error(
        "This participant attempted to modify the Mesh \"{}\" while locked. "
        "Mesh modification is only allowed before calling initialize().",
        meshNameStr));
  }
  const int expectedDims = _impl->configData.meshes.at(meshNameStr).dimensions;
  // Expect coordinates as expectedDims * N entries for N VertexIDs
  if (coordinates.size() != ids.size() * expectedDims) {
    throw precice::Error(precice::utils::format_or_error("setMeshVertices() was called with {} vertices and {} coordinates ({}D), but needs {} coordinates ({} x {}).", ids.size(), coordinates.size(), expectedDims, ids.size() * expectedDims, ids.size(), expectedDims));
  }
  auto &count = _impl->runtimeState.meshVertexCounts[meshNameStr];
  for (std::size_t i = 0; i < ids.size(); ++i) {
    ids[i] = VertexID(count + i);
  }
  count += ids.size();
}

void Participant::setMeshEdge(
    ::precice::string_view meshName,
    VertexID               first,
    VertexID               second)
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  std::string                           meshNameStr = toStr(meshName);
  _impl->requireMeshUsed(meshNameStr);
  if (_impl->runtimeState.lockedMeshes.count(meshNameStr)) {
    throw precice::Error(precice::utils::format_or_error(
        "This participant attempted to modify the Mesh \"{}\" while locked. "
        "Mesh modification is only allowed before calling initialize().",
        meshNameStr));
  }
  if (!_impl->configData.meshes.at(meshNameStr).requiresConnectivity) {
    return;
  }
  auto countIt = _impl->runtimeState.meshVertexCounts.find(meshNameStr);
  auto count   = (countIt != _impl->runtimeState.meshVertexCounts.end()) ? static_cast<VertexID>(countIt->second) : VertexID(0);
  if (first < 0 || first >= count) {
    throw precice::Error(precice::utils::format_or_error("setMeshEdge() was called with an invalid vertex ID {}.", first));
  }
  if (second < 0 || second >= count) {
    throw precice::Error(precice::utils::format_or_error("setMeshEdge() was called with an invalid vertex ID {}.", second));
  }
  // no-op
}

void Participant::setMeshEdges(
    ::precice::string_view        meshName,
    precice::span<const VertexID> ids)
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  std::string                           meshNameStr = toStr(meshName);
  _impl->requireMeshUsed(meshNameStr);
  if (_impl->runtimeState.lockedMeshes.count(meshNameStr)) {
    throw precice::Error(precice::utils::format_or_error(
        "This participant attempted to modify the Mesh \"{}\" while locked. "
        "Mesh modification is only allowed before calling initialize().",
        meshNameStr));
  }
  // Like the real library, all further validation is skipped when the mesh does not
  // require connectivity (the call is a no-op then).
  if (!_impl->configData.meshes.at(meshNameStr).requiresConnectivity) {
    return;
  }
  if (ids.size() % 2 != 0) {
    throw precice::Error(precice::utils::format_or_error("setMeshEdges() was called with {} vertices, but needs an even number of vertices (2 per edge).", ids.size()));
  }
  auto countIt = _impl->runtimeState.meshVertexCounts.find(meshNameStr);
  auto count   = (countIt != _impl->runtimeState.meshVertexCounts.end()) ? static_cast<VertexID>(countIt->second) : VertexID(0);
  for (std::size_t i = 0; i < ids.size(); ++i) {
    if (ids[i] < 0 || ids[i] >= count) {
      throw precice::Error(precice::utils::format_or_error("setMeshEdges() was called with an invalid vertex ID {} at index {}.", ids[i], i));
    }
  }
  // no-op
}

void Participant::setMeshTriangle(
    ::precice::string_view meshName,
    VertexID               first,
    VertexID               second,
    VertexID               third)
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  std::string                           meshNameStr = toStr(meshName);
  _impl->requireMeshUsed(meshNameStr);
  if (_impl->runtimeState.lockedMeshes.count(meshNameStr)) {
    throw precice::Error(precice::utils::format_or_error(
        "This participant attempted to modify the Mesh \"{}\" while locked. "
        "Mesh modification is only allowed before calling initialize().",
        meshNameStr));
  }
  if (!_impl->configData.meshes.at(meshNameStr).requiresConnectivity) {
    return;
  }
  auto countIt = _impl->runtimeState.meshVertexCounts.find(meshNameStr);
  auto count   = (countIt != _impl->runtimeState.meshVertexCounts.end()) ? static_cast<VertexID>(countIt->second) : VertexID(0);
  if (first < 0 || first >= count) {
    throw precice::Error(precice::utils::format_or_error("setMeshTriangle() was called with an invalid vertex ID {}.", first));
  }
  if (second < 0 || second >= count) {
    throw precice::Error(precice::utils::format_or_error("setMeshTriangle() was called with an invalid vertex ID {}.", second));
  }
  if (third < 0 || third >= count) {
    throw precice::Error(precice::utils::format_or_error("setMeshTriangle() was called with an invalid vertex ID {}.", third));
  }
  if (first == second || first == third || second == third) {
    throw precice::Error(precice::utils::format_or_error("setMeshTriangle() was called with repeated Vertex IDs ({}, {}, {}).", first, second, third));
  }
  // no-op
}

void Participant::setMeshTriangles(
    ::precice::string_view        meshName,
    precice::span<const VertexID> ids)
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  std::string                           meshNameStr = toStr(meshName);
  _impl->requireMeshUsed(meshNameStr);
  if (_impl->runtimeState.lockedMeshes.count(meshNameStr)) {
    throw precice::Error(precice::utils::format_or_error(
        "This participant attempted to modify the Mesh \"{}\" while locked. "
        "Mesh modification is only allowed before calling initialize().",
        meshNameStr));
  }
  if (!_impl->configData.meshes.at(meshNameStr).requiresConnectivity) {
    return;
  }
  if (ids.size() % 3 != 0) {
    throw precice::Error(precice::utils::format_or_error("setMeshTriangles() was called with {} vertices, but needs a number of vertices divisible by 3 (3 per triangle).", ids.size()));
  }
  auto countIt = _impl->runtimeState.meshVertexCounts.find(meshNameStr);
  auto count   = (countIt != _impl->runtimeState.meshVertexCounts.end()) ? static_cast<VertexID>(countIt->second) : VertexID(0);
  for (std::size_t i = 0; i < ids.size(); ++i) {
    if (ids[i] < 0 || ids[i] >= count) {
      throw precice::Error(precice::utils::format_or_error("setMeshTriangles() was called with an invalid vertex ID {} at index {}.", ids[i], i));
    }
  }
  // no-op
}

void Participant::setMeshQuad(
    ::precice::string_view meshName,
    VertexID               first,
    VertexID               second,
    VertexID               third,
    VertexID               fourth)
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  std::string                           meshNameStr = toStr(meshName);
  _impl->requireMeshUsed(meshNameStr);
  if (_impl->runtimeState.lockedMeshes.count(meshNameStr)) {
    throw precice::Error(precice::utils::format_or_error(
        "This participant attempted to modify the Mesh \"{}\" while locked. "
        "Mesh modification is only allowed before calling initialize().",
        meshNameStr));
  }
  if (!_impl->configData.meshes.at(meshNameStr).requiresConnectivity) {
    return;
  }
  auto countIt = _impl->runtimeState.meshVertexCounts.find(meshNameStr);
  auto count   = (countIt != _impl->runtimeState.meshVertexCounts.end()) ? static_cast<VertexID>(countIt->second) : VertexID(0);
  if (first < 0 || first >= count) {
    throw precice::Error(precice::utils::format_or_error("setMeshQuad() was called with an invalid vertex ID {}.", first));
  }
  if (second < 0 || second >= count) {
    throw precice::Error(precice::utils::format_or_error("setMeshQuad() was called with an invalid vertex ID {}.", second));
  }
  if (third < 0 || third >= count) {
    throw precice::Error(precice::utils::format_or_error("setMeshQuad() was called with an invalid vertex ID {}.", third));
  }
  if (fourth < 0 || fourth >= count) {
    throw precice::Error(precice::utils::format_or_error("setMeshQuad() was called with an invalid vertex ID {}.", fourth));
  }
  if (first == second || first == third || first == fourth || second == third || second == fourth || third == fourth) {
    throw precice::Error("The four vertex ID's are not unique. Please check that the vertices that form the quad are correct.");
  }
  // no-op
}

void Participant::setMeshQuads(
    ::precice::string_view        meshName,
    precice::span<const VertexID> ids)
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  std::string                           meshNameStr = toStr(meshName);
  _impl->requireMeshUsed(meshNameStr);
  if (_impl->runtimeState.lockedMeshes.count(meshNameStr)) {
    throw precice::Error(precice::utils::format_or_error(
        "This participant attempted to modify the Mesh \"{}\" while locked. "
        "Mesh modification is only allowed before calling initialize().",
        meshNameStr));
  }
  if (!_impl->configData.meshes.at(meshNameStr).requiresConnectivity) {
    return;
  }
  if (ids.size() % 4 != 0) {
    throw precice::Error(precice::utils::format_or_error("setMeshQuads() was called with {} vertices, but needs a number of vertices divisible by 4 (4 per quad).", ids.size()));
  }
  auto countIt = _impl->runtimeState.meshVertexCounts.find(meshNameStr);
  auto count   = (countIt != _impl->runtimeState.meshVertexCounts.end()) ? static_cast<VertexID>(countIt->second) : VertexID(0);
  for (std::size_t i = 0; i < ids.size(); ++i) {
    if (ids[i] < 0 || ids[i] >= count) {
      throw precice::Error(precice::utils::format_or_error("setMeshQuads() was called with an invalid vertex ID {} at index {}.", ids[i], i));
    }
  }
  for (std::size_t i = 0; i < ids.size() / 4; ++i) {
    auto a = ids[4 * i], b = ids[4 * i + 1], c = ids[4 * i + 2], d = ids[4 * i + 3];
    if (a == b || a == c || a == d || b == c || b == d || c == d) {
      throw precice::Error(precice::utils::format_or_error("The four vertex ID's of the quad nr {} are not unique. Please check that the vertices that form the quad are correct.", i));
    }
  }
  // no-op
}

void Participant::setMeshTetrahedron(
    ::precice::string_view meshName,
    VertexID               first,
    VertexID               second,
    VertexID               third,
    VertexID               fourth)
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  std::string                           meshNameStr = toStr(meshName);
  _impl->requireMeshUsed(meshNameStr);
  if (_impl->runtimeState.lockedMeshes.count(meshNameStr)) {
    throw precice::Error(precice::utils::format_or_error(
        "This participant attempted to modify the Mesh \"{}\" while locked. "
        "Mesh modification is only allowed before calling initialize().",
        meshNameStr));
  }
  const auto &meshInfo = _impl->configData.meshes.at(meshNameStr);
  if (meshInfo.dimensions != 3) {
    throw precice::Error("setMeshTetrahedron is only possible for 3D meshes. "
                         "Please set the mesh dimension to 3 in the preCICE configuration file.");
  }
  if (!meshInfo.requiresConnectivity) {
    return;
  }
  auto countIt = _impl->runtimeState.meshVertexCounts.find(meshNameStr);
  auto count   = (countIt != _impl->runtimeState.meshVertexCounts.end()) ? static_cast<VertexID>(countIt->second) : VertexID(0);
  if (first < 0 || first >= count) {
    throw precice::Error(precice::utils::format_or_error("setMeshTetrahedron() was called with an invalid vertex ID {}.", first));
  }
  if (second < 0 || second >= count) {
    throw precice::Error(precice::utils::format_or_error("setMeshTetrahedron() was called with an invalid vertex ID {}.", second));
  }
  if (third < 0 || third >= count) {
    throw precice::Error(precice::utils::format_or_error("setMeshTetrahedron() was called with an invalid vertex ID {}.", third));
  }
  if (fourth < 0 || fourth >= count) {
    throw precice::Error(precice::utils::format_or_error("setMeshTetrahedron() was called with an invalid vertex ID {}.", fourth));
  }
  if (first == second || first == third || first == fourth || second == third || second == fourth || third == fourth) {
    throw precice::Error(precice::utils::format_or_error("setMeshTetrahedron() was called with vertices [{}, {}, {}, {}], but some vertices are equal.", first, second, third, fourth));
  }
  // no-op
}

void Participant::setMeshTetrahedra(
    ::precice::string_view        meshName,
    precice::span<const VertexID> ids)
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  std::string                           meshNameStr = toStr(meshName);
  _impl->requireMeshUsed(meshNameStr);
  if (_impl->runtimeState.lockedMeshes.count(meshNameStr)) {
    throw precice::Error(precice::utils::format_or_error(
        "This participant attempted to modify the Mesh \"{}\" while locked. "
        "Mesh modification is only allowed before calling initialize().",
        meshNameStr));
  }
  const auto &meshInfo = _impl->configData.meshes.at(meshNameStr);
  if (meshInfo.dimensions != 3) {
    throw precice::Error("setMeshTetrahedron is only possible for 3D meshes. "
                         "Please set the mesh dimension to 3 in the preCICE configuration file.");
  }
  if (!meshInfo.requiresConnectivity) {
    return;
  }
  if (ids.size() % 4 != 0) {
    throw precice::Error(precice::utils::format_or_error("setMeshTetrahedra() was called with {} vertices, but needs a number of vertices divisible by 4 (4 per tetrahedron).", ids.size()));
  }
  auto countIt = _impl->runtimeState.meshVertexCounts.find(meshNameStr);
  auto count   = (countIt != _impl->runtimeState.meshVertexCounts.end()) ? static_cast<VertexID>(countIt->second) : VertexID(0);
  for (std::size_t i = 0; i < ids.size(); ++i) {
    if (ids[i] < 0 || ids[i] >= count) {
      throw precice::Error(precice::utils::format_or_error("setMeshTetrahedra() was called with an invalid vertex ID {} at index {}.", ids[i], i));
    }
  }
  // no-op
}

// Data access

bool Participant::requiresInitialData()
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (_impl->runtimeState.initialized) {
    throw precice::Error("requiresInitialData() has to be called before initialize().");
  }

  if (_impl->runtimeState.initialDataRequired) {
    _impl->runtimeState.initialDataFulfilled = true;
    return true;
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
  if (_impl->runtimeState.finalized) {
    throw precice::Error("writeData() cannot be called after finalize().");
  }
  if (_impl->runtimeState.initialized && !_impl->runtimeState.couplingOngoing) {
    throw precice::Error("Calling writeData(...) is forbidden if coupling is not ongoing, because the data you are trying to write will not be used anymore. "
                         "You can fix this by always calling writeData(...) before the advance(...) call in your simulation loop "
                         "or by using Participant::isCouplingOngoing() to implement a safeguard.");
  }

  std::string meshNameStr  = toStr(meshName);
  std::string dataNameStr  = toStr(dataName);
  int         expectedDims = 1;

  // Validate against config if parsed
  if (_impl->runtimeState.configParsed) {
    _impl->requireMeshUsed(meshNameStr);
    bool found = false;

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

  if (ids.empty() && values.empty()) {
    return;
  }

  // Validate vertex IDs
  auto countIt = _impl->runtimeState.meshVertexCounts.find(meshNameStr);
  auto count   = (countIt != _impl->runtimeState.meshVertexCounts.end()) ? static_cast<VertexID>(countIt->second) : VertexID(0);
  for (std::size_t i = 0; i < ids.size(); ++i) {
    if (ids[i] < 0 || ids[i] >= count) {
      throw precice::Error(precice::utils::format_or_error(
          "Cannot write data \"{}\" to mesh \"{}\" due to invalid Vertex ID at vertices[{}]. "
          "Please make sure you only use the results from calls to setMeshVertex/Vertices().",
          dataNameStr, meshNameStr, i));
    }
  }

  // Validate data values (no NaN or Inf). Real preCICE only performs this check in debug
  // builds (PRECICE_VALIDATE_DATA is a no-op under NDEBUG), so we mirror that to avoid being
  // stricter than the release library that tutorials run against.
#ifndef NDEBUG
  if (!std::all_of(values.begin(), values.end(), [](double v) { return std::isfinite(v); })) {
    throw precice::Error("One of the given data values is either plus or minus infinity or NaN.");
  }
#endif

  // Store written data per mesh:data pair; also keep the globally last-written
  // data as a fallback for reads of data that is never written by this participant.
  _impl->runtimeState.writeBuffers[meshNameStr + ":" + dataNameStr].assign(values.begin(), values.end());
  _impl->runtimeState.lastWriteBuffer.assign(values.begin(), values.end());
  // Track last write size for this mesh+data to validate future reads
  _impl->runtimeState.lastWriteSizes[meshNameStr + ":" + dataNameStr] = values.size();
}

void Participant::readData(
    ::precice::string_view        meshName,
    ::precice::string_view        dataName,
    precice::span<const VertexID> ids,
    double                        relativeReadTime,
    precice::span<double>         values) const
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  const std::size_t                     n = values.size();
  if (_impl->runtimeState.finalized) {
    throw precice::Error("readData(...) cannot be called after finalize().");
  }
  if (!_impl->runtimeState.initialized) {
    throw precice::Error("readData(...) cannot be called before initialize().");
  }
  if (relativeReadTime < 0.0) {
    throw precice::Error("readData(...) cannot sample data before the current time.");
  }
  // Maximum sample time = time until end of window, truncated by max-time, or 0 once coupling
  // ended (mirrors getNextTimeStepMaxSize). Compared with tolerance like math::smallerEquals.
  const double maxDt = _impl->nextTimeStepMaxSize();
  if (relativeReadTime - maxDt > timeTolerance) {
    throw precice::Error("readData(...) cannot sample data outside of current time window.");
  }
  if (!_impl->runtimeState.couplingOngoing && relativeReadTime != 0.0) {
    throw precice::Error(precice::utils::format_or_error(
        "Calling readData(...) with relativeReadTime = {} is forbidden if coupling is not ongoing. "
        "If coupling finished, only data for relativeReadTime = 0 is available. "
        "Please always use precice.getMaxTimeStepSize() to obtain the maximum allowed relativeReadTime.",
        relativeReadTime));
  }

  std::string meshNameStr  = toStr(meshName);
  std::string dataNameStr  = toStr(dataName);
  int         expectedDims = 1;

  if (_impl->runtimeState.configParsed) {
    _impl->requireMeshUsed(meshNameStr);
  }

  // Check mesh is not in reset state
  if (_impl->runtimeState.initialized && _impl->runtimeState.lockedMeshes.count(meshNameStr) == 0) {
    throw precice::Error(precice::utils::format_or_error(
        "Cannot read from mesh \"{}\" after it has been reset. Please read data before calling resetMesh().",
        meshNameStr));
  }

  // Validate against config if parsed
  if (_impl->runtimeState.configParsed) {
    bool found = false;

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

  if (ids.empty() && values.empty()) {
    return;
  }

  // Validate vertex IDs
  {
    auto countIt = _impl->runtimeState.meshVertexCounts.find(meshNameStr);
    auto count   = (countIt != _impl->runtimeState.meshVertexCounts.end()) ? static_cast<VertexID>(countIt->second) : VertexID(0);
    for (std::size_t i = 0; i < ids.size(); ++i) {
      if (ids[i] < 0 || ids[i] >= count) {
        throw precice::Error(precice::utils::format_or_error(
            "Cannot read data \"{}\" from mesh \"{}\" due to invalid Vertex ID at vertices[{}]. "
            "Please make sure you only use the results from calls to setMeshVertex/Vertices().",
            dataNameStr, meshNameStr, i));
      }
    }
  }

  // Determine data mode from mock config
  std::string                     key        = meshNameStr + ":" + dataNameStr;
  impl::ParticipantImpl::DataMode mode       = _impl->mockConfig.defaultMode;
  double                          scalarMult = _impl->mockConfig.defaultScalarMultiplier;
  double                          randLower  = _impl->mockConfig.defaultRandomLower;
  double                          randUpper  = _impl->mockConfig.defaultRandomUpper;
  uint32_t                        randSeed   = _impl->mockConfig.defaultRandomSeed;
  std::vector<double>             vectorMult = _impl->mockConfig.defaultVectorMultiplier;

  auto mockIt = _impl->mockConfig.dataConfigs.find(key);
  if (mockIt != _impl->mockConfig.dataConfigs.end()) {
    mode       = mockIt->second.mode;
    scalarMult = mockIt->second.scalarMultiplier;
    randLower  = mockIt->second.randomLower;
    randUpper  = mockIt->second.randomUpper;
    randSeed   = mockIt->second.randomSeed;
    vectorMult = mockIt->second.vectorMultiplier;
  }

  // Validate read size against last write for the same mesh/data (if available)
  auto sizeIt = _impl->runtimeState.lastWriteSizes.find(meshNameStr + ":" + dataNameStr);
  if (sizeIt != _impl->runtimeState.lastWriteSizes.end() && sizeIt->second != n) {
    throw precice::Error(precice::utils::format_or_error(
        "readData() was called with {} values for data '{}' on mesh '{}', but the last write had {} values for this mesh/data.",
        n, dataNameStr, meshNameStr, sizeIt->second));
  }

  // Resolve the buffer to echo back: exact mesh:data buffer first, then the same
  // data name written on another mesh, then the globally last-written data.
  const std::vector<double> *sourceBuffer = nullptr;
  {
    auto bufIt = _impl->runtimeState.writeBuffers.find(key);
    if (bufIt != _impl->runtimeState.writeBuffers.end() && !bufIt->second.empty()) {
      sourceBuffer = &bufIt->second;
    } else {
      const std::string suffix = ":" + dataNameStr;
      for (const auto &entry : _impl->runtimeState.writeBuffers) {
        if (entry.first.size() > suffix.size() &&
            entry.first.compare(entry.first.size() - suffix.size(), suffix.size(), suffix) == 0 &&
            !entry.second.empty()) {
          sourceBuffer = &entry.second;
          break;
        }
      }
      if (sourceBuffer == nullptr && !_impl->runtimeState.lastWriteBuffer.empty()) {
        sourceBuffer = &_impl->runtimeState.lastWriteBuffer;
      }
    }
  }

  // Apply the selected mode
  switch (mode) {
  case impl::ParticipantImpl::DataMode::Random: {
    // Mode 1: Random data with configurable bounds and optional seed
    if (!std::isfinite(randLower) || !std::isfinite(randUpper) || randUpper <= randLower) {
      throw precice::Error(precice::utils::format_or_error(
          "Invalid random bounds for data '{}' on mesh '{}': upper={} must be greater than lower={} and both must be finite.",
          dataNameStr, meshNameStr, randUpper, randLower));
    }
    // Use custom seed if provided (non-zero), otherwise use default seed
    uint32_t                               effectiveSeed = (randSeed != 0) ? randSeed : (static_cast<uint32_t>(_impl->runtimeState.seed + static_cast<uint32_t>(_impl->runtimeState.currentStep)));
    std::mt19937                           gen(effectiveSeed);
    std::uniform_real_distribution<double> dist(randLower, randUpper);
    for (std::size_t i = 0; i < n; ++i) {
      values[i] = dist(gen);
    }
    break;
  }

  case impl::ParticipantImpl::DataMode::Buffer: {
    // Mode 2: Return buffered write data (cyclic if shorter)
    if (sourceBuffer != nullptr) {
      const auto &buffer = *sourceBuffer;
      for (std::size_t i = 0; i < n; ++i) {
        values[i] = buffer[i % buffer.size()];
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
    // Mode 3: Return buffered write data with strict scaling (one multiplier per component/value)
    if (sourceBuffer != nullptr) {
      const auto &buffer = *sourceBuffer;
      if (!vectorMult.empty()) {
        for (std::size_t i = 0; i < n; ++i) {
          values[i] = buffer[i % buffer.size()] * vectorMult[i % vectorMult.size()];
        }
      } else {
        // Scalar multiplication
        for (std::size_t i = 0; i < n; ++i) {
          values[i] = buffer[i % buffer.size()] * scalarMult;
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
    ::precice::string_view      meshName,
    ::precice::string_view      dataName,
    precice::span<const double> coordinates,
    precice::span<const double> values)
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (!_impl->configData.allowsExperimental) {
    throw precice::Error("You called the API function \"writeAndMapData\", which is part of the experimental API. "
                         "You may unlock the full API by specifying <precice-configuration experimental=\"true\" ... > in the configuration. Please be aware that experimental features may change in any future version (even minor or bugfix).");
  }
  if (_impl->runtimeState.finalized) {
    throw precice::Error("writeAndMapData() cannot be called after finalize().");
  }
  if (!_impl->runtimeState.initialized) {
    throw precice::Error("writeAndMapData() cannot be called before initialize(), because the mesh to map onto hasn't been received yet.");
  }
  if (!_impl->runtimeState.couplingOngoing) {
    throw precice::Error("writeAndMapData() is forbidden if coupling is not ongoing.");
  }
  std::string meshNameStr = toStr(meshName);
  std::string dataNameStr = toStr(dataName);
  // Check received mesh with api-access
  auto meshIt = _impl->configData.meshes.find(meshNameStr);
  if (meshIt == _impl->configData.meshes.end() || !meshIt->second.received || !meshIt->second.apiAccess) {
    throw precice::Error(precice::utils::format_or_error(
        "This participant attempted to map and write data (via \"writeAndMapData\") to mesh \"{0}\", "
        "but mesh \"{0}\" is either not a received mesh or its api access was not enabled in the configuration. "
        "writeAndMapData({0}, ...) is only valid for (<receive-mesh name=\"{0}\" ... api-access=\"true\"/>).",
        meshNameStr));
  }
  // Check JIT mapping
  {
    bool hasJit = false;
    for (const auto &m : _impl->configData.mappings) {
      if (m.direction == "write" && m.isJit && m.toMesh == meshNameStr) {
        hasJit = true;
        break;
      }
    }
    if (!hasJit) {
      throw precice::Error(precice::utils::format_or_error(
          "The function \"writeAndMapData\" was called on mesh \"{0}\", but no matching just-in-time mapping was configured. "
          "Please define a mapping in write direction to the mesh \"{0}\" and omit the \"from\" attribute from the definition.",
          meshNameStr));
    }
  }
  if (_impl->runtimeState.meshAccessRegions.find(meshNameStr) == _impl->runtimeState.meshAccessRegions.end()) {
    throw precice::Error(precice::utils::format_or_error(
        "The function \"writeAndMapData\" was called on mesh \"{}\", but no access region was defined. Please define an access region using \"setMeshAccessRegion()\" before calling \"writeAndMapData()\".",
        meshNameStr));
  }
  const int dim = getMeshDimensions(meshName);
  if (dim <= 0) {
    throw precice::Error(precice::utils::format_or_error("writeAndMapData() was called for unknown mesh '{}'.", meshNameStr));
  }
  if (coordinates.size() % static_cast<std::size_t>(dim) != 0) {
    throw precice::Error(precice::utils::format_or_error("writeAndMapData() was called with {} coordinates, but needs a multiple of {} coordinates.", coordinates.size(), dim));
  }
  const std::size_t nVertices      = coordinates.size() / static_cast<std::size_t>(dim);
  const int         dataDims       = getDataDimensions(meshName, dataName);
  const std::size_t expectedValues = nVertices * static_cast<std::size_t>(dataDims);
  if (values.size() != expectedValues) {
    throw precice::Error(precice::utils::format_or_error(
        "Input sizes are inconsistent attempting to write {}D data to mesh. "
        "You passed {} vertex indices and {} data components, but we expected {} data components ({} x {}).",
        dataDims, nVertices, values.size(), expectedValues, dataDims, nVertices));
  }
  // no-op
}

void Participant::mapAndReadData(
    ::precice::string_view      meshName,
    ::precice::string_view      dataName,
    precice::span<const double> coordinates,
    double                      relativeReadTime,
    precice::span<double>       values) const
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (!_impl->configData.allowsExperimental) {
    throw precice::Error("You called the API function \"mapAndReadData\", which is part of the experimental API. "
                         "You may unlock the full API by specifying <precice-configuration experimental=\"true\" ... > in the configuration. Please be aware that experimental features may change in any future version (even minor or bugfix).");
  }
  if (_impl->runtimeState.finalized) {
    throw precice::Error("mapAndReadData() cannot be called after finalize().");
  }
  if (!_impl->runtimeState.initialized) {
    throw precice::Error("initialize() has to be called before mapAndReadData().");
  }
  if (!_impl->runtimeState.couplingOngoing && relativeReadTime != 0.0) {
    throw precice::Error(precice::utils::format_or_error(
        "Calling mapAndReadData() with relativeReadTime = {} is forbidden if coupling is not ongoing. If coupling finished, only data for relativeReadTime = 0 is available.",
        relativeReadTime));
  }
  std::string meshNameStr = toStr(meshName);
  // Check received mesh with api-access
  auto meshIt = _impl->configData.meshes.find(meshNameStr);
  if (meshIt == _impl->configData.meshes.end() || !meshIt->second.received || !meshIt->second.apiAccess) {
    throw precice::Error(precice::utils::format_or_error(
        "This participant attempted to map and read data (via \"mapAndReadData\") from mesh \"{0}\", "
        "but mesh \"{0}\" is either not a received mesh or its api access was not enabled in the configuration. "
        "mapAndReadData({0}, ...) is only valid for (<receive-mesh name=\"{0}\" ... api-access=\"true\"/>).",
        meshNameStr));
  }
  // Check JIT mapping
  {
    bool hasJit = false;
    for (const auto &m : _impl->configData.mappings) {
      if (m.direction == "read" && m.isJit && m.fromMesh == meshNameStr) {
        hasJit = true;
        break;
      }
    }
    if (!hasJit) {
      throw precice::Error(precice::utils::format_or_error(
          "The function \"mapAndReadData\" was called on mesh \"{0}\", but no matching just-in-time mapping was configured. "
          "Please define a mapping in read direction from the mesh \"{0}\" and omit the \"to\" attribute from the definition.",
          meshNameStr));
    }
  }
  if (_impl->runtimeState.meshAccessRegions.find(meshNameStr) == _impl->runtimeState.meshAccessRegions.end()) {
    throw precice::Error(precice::utils::format_or_error(
        "The function \"mapAndReadData\" was called on mesh \"{}\", but no access region was defined. Please define an access region using \"setMeshAccessRegion()\" before calling \"mapAndReadData()\".",
        meshNameStr));
  }
  const int dim = getMeshDimensions(meshName);
  if (dim <= 0) {
    throw precice::Error(precice::utils::format_or_error("mapAndReadData() was called for unknown mesh '{}'.", meshNameStr));
  }
  if (coordinates.size() % static_cast<std::size_t>(dim) != 0) {
    throw precice::Error(precice::utils::format_or_error("mapAndReadData() was called with {} coordinates, but needs a multiple of {} coordinates.", coordinates.size(), dim));
  }
  const std::size_t nVertices      = coordinates.size() / static_cast<std::size_t>(dim);
  const int         dataDims       = getDataDimensions(meshName, dataName);
  const std::size_t expectedValues = nVertices * static_cast<std::size_t>(dataDims);
  if (values.size() != expectedValues) {
    throw precice::Error(precice::utils::format_or_error(
        "Input/Output sizes are inconsistent attempting to read {}D data from mesh. "
        "You passed {} vertex indices and {} data components, but we expected {} data components ({} x {}).",
        dataDims, nVertices, values.size(), expectedValues, dataDims, nVertices));
  }
  // no-op
}

// Direct access

void Participant::setMeshAccessRegion(
    ::precice::string_view      meshName,
    precice::span<const double> boundingBox) const
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  std::string                           meshNameStr = toStr(meshName);
  auto                                  meshIt      = _impl->configData.meshes.find(meshNameStr);
  if (meshIt == _impl->configData.meshes.end() || !meshIt->second.received || !meshIt->second.apiAccess) {
    throw precice::Error(precice::utils::format_or_error(
        "This participant attempted to set an access region (via \"setMeshAccessRegion\") on mesh \"{0}\", "
        "but mesh \"{0}\" is either not a received mesh or its api access was not enabled in the configuration. "
        "setMeshAccessRegion(...) is only valid for (<receive-mesh name=\"{0}\" ... api-access=\"true\"/>).",
        meshNameStr));
  }
  if (_impl->runtimeState.finalized) {
    throw precice::Error("setMeshAccessRegion() cannot be called after finalize().");
  }
  if (_impl->runtimeState.initialized) {
    throw precice::Error("setMeshAccessRegion() needs to be called before initialize().");
  }
  if (_impl->runtimeState.meshAccessRegions.find(meshNameStr) != _impl->runtimeState.meshAccessRegions.end()) {
    throw precice::Error(precice::utils::format_or_error(
        "A mesh access region was already defined for mesh '{}'. setMeshAccessRegion may only be called once per mesh.",
        meshNameStr));
  }
  const auto meshDimensions = getMeshDimensions(meshName);
  if (meshDimensions <= 0) {
    throw precice::Error(precice::utils::format_or_error(
        "setMeshAccessRegion() was called for unknown mesh '{}'.", meshNameStr));
  }
  const auto expectedSize = static_cast<std::size_t>(2 * meshDimensions);
  if (boundingBox.size() != expectedSize) {
    throw precice::Error(precice::utils::format_or_error(
        "setMeshAccessRegion() was called with {} bounding box coordinates, but expects {} coordinates (min and max for each dimension in {}D).",
        boundingBox.size(), expectedSize, meshDimensions));
  }
  for (int d = 0; d < meshDimensions; ++d) {
    if (boundingBox[2 * d] > boundingBox[2 * d + 1]) {
      throw precice::Error("Your bounding box is ill defined, i.e. it has a negative volume. The required format is [x_min, x_max...]");
    }
  }
  _impl->runtimeState.meshAccessRegions[meshNameStr] = std::vector<double>(boundingBox.begin(), boundingBox.end());

  // The partner's mesh is not available in the mock, so synthesize a small set of
  // vertices inside the access region. Access regions are typically a thin box that
  // straddles the coupling interface: wide along the interface, thin along its normal.
  // Spreading points along the box diagonal pushes them across the thin dimension and
  // past the interface, landing outside the partner's domain — so solvers that locate
  // these points on their own mesh (e.g. nutils) fail (LocateError). Instead, spread the
  // points only along the widest box dimension and place every other dimension at the box
  // center: the thin dimension then lands on the interface plane, where the partner's
  // mesh actually lives, keeping the points locatable.
  constexpr std::size_t kSyntheticVertices = 10;
  // Keep the spread within the interior band of the box: the access region is usually a
  // superset of the partner's domain, so points near the ends often fall just outside it.
  constexpr double kMargin = 0.25; // map t into [kMargin, 1 - kMargin]

  int    widestDim   = 0;
  double widestWidth = -1.0;
  for (int d = 0; d < meshDimensions; ++d) {
    const double width = boundingBox[2 * d + 1] - boundingBox[2 * d];
    if (width > widestWidth) {
      widestWidth = width;
      widestDim   = d;
    }
  }

  std::vector<double> coords;
  coords.reserve(kSyntheticVertices * meshDimensions);
  for (std::size_t v = 0; v < kSyntheticVertices; ++v) {
    const double raw = (kSyntheticVertices == 1) ? 0.5 : static_cast<double>(v) / (kSyntheticVertices - 1);
    const double t   = kMargin + raw * (1.0 - 2.0 * kMargin);
    for (int d = 0; d < meshDimensions; ++d) {
      const double lo   = boundingBox[2 * d];
      const double hi   = boundingBox[2 * d + 1];
      const double frac = (d == widestDim) ? t : 0.5; // spread along widest, center the rest
      coords.push_back(lo + frac * (hi - lo));
    }
  }
  _impl->runtimeState.directAccessCoords[meshNameStr] = std::move(coords);
  _impl->runtimeState.meshVertexCounts[meshNameStr]   = kSyntheticVertices;
}

void Participant::getMeshVertexIDsAndCoordinates(
    ::precice::string_view  meshName,
    precice::span<VertexID> ids,
    precice::span<double>   coordinates) const
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  std::string                           meshNameStr = toStr(meshName);
  auto                                  meshIt      = _impl->configData.meshes.find(meshNameStr);
  if (meshIt == _impl->configData.meshes.end() || !meshIt->second.received || !meshIt->second.apiAccess) {
    throw precice::Error(precice::utils::format_or_error(
        "This participant attempted to get mesh vertex IDs and coordinates (via \"getMeshVertexIDsAndCoordinates\") from mesh \"{0}\", "
        "but mesh \"{0}\" is either not a received mesh or its api access was not enabled in the configuration. "
        "getMeshVertexIDsAndCoordinates(...) is only valid for (<receive-mesh name=\"{0}\" ... api-access=\"true\"/>).",
        meshNameStr));
  }
  if (!_impl->runtimeState.initialized) {
    throw precice::Error("initialize() has to be called before getMeshVertexIDsAndCoordinates().");
  }
  if (_impl->runtimeState.meshAccessRegions.find(meshNameStr) == _impl->runtimeState.meshAccessRegions.end()) {
    throw precice::Error("Please define an access region using \"setMeshAccessRegion()\" before calling \"getMeshVertexIDsAndCoordinates()\".");
  }
  const auto meshDimensions = getMeshDimensions(meshName);
  if (meshDimensions <= 0) {
    throw precice::Error(precice::utils::format_or_error(
        "getMeshVertexIDsAndCoordinates() was called for unknown mesh '{}'.", meshNameStr));
  }
  const auto expectedCoordinates = ids.size() * static_cast<std::size_t>(meshDimensions);
  if (coordinates.size() != expectedCoordinates) {
    throw precice::Error(precice::utils::format_or_error("getMeshVertexIDsAndCoordinates() was called with {} vertices and {} coordinates, but needs {} coordinates ({} x {}).", ids.size(), coordinates.size(), expectedCoordinates, ids.size(), meshDimensions));
  }

  const auto vertexCount = static_cast<std::size_t>(getMeshVertexSize(meshName));
  if (ids.size() != vertexCount) {
    throw precice::Error(precice::utils::format_or_error("getMeshVertexIDsAndCoordinates() was called with {} vertex IDs, but the mesh contains {} vertices.", ids.size(), vertexCount));
  }

  for (std::size_t i = 0; i < ids.size(); ++i) {
    ids[i] = static_cast<VertexID>(i);
  }
  auto coordIt = _impl->runtimeState.directAccessCoords.find(meshNameStr);
  if (coordIt != _impl->runtimeState.directAccessCoords.end() && coordIt->second.size() == coordinates.size()) {
    std::copy(coordIt->second.begin(), coordIt->second.end(), coordinates.begin());
  } else {
    std::fill(coordinates.begin(), coordinates.end(), 0.0);
  }
}

// Gradient data (experimental)

bool Participant::requiresGradientDataFor(::precice::string_view meshName,
                                          ::precice::string_view dataName) const
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  std::string                           meshNameStr = toStr(meshName);
  std::string                           dataNameStr = toStr(dataName);

  _impl->requireMeshUsed(meshNameStr);
  bool dataKnown = false;

  // Read data never requires gradients
  for (const auto &dataInfo : _impl->configData.dataItems) {
    if (dataInfo.name == dataNameStr && dataInfo.meshName == meshNameStr) {
      dataKnown = true;
      if (dataInfo.isRead) {
        return false;
      }
    }
  }
  if (!dataKnown) {
    throw precice::Error(precice::utils::format_or_error(
        "The given data \"{}\" is unknown to preCICE mesh \"{}\".", dataNameStr, meshNameStr));
  }

  // Check if any write mapping from this mesh requires gradient data
  for (const auto &m : _impl->configData.mappings) {
    if (m.direction == "write" && m.fromMesh == meshNameStr && m.requiresGradient) {
      return true;
    }
  }
  return false;
}

void Participant::writeGradientData(
    ::precice::string_view        meshName,
    ::precice::string_view        dataName,
    precice::span<const VertexID> ids,
    precice::span<const double>   gradients)
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (!_impl->configData.allowsExperimental) {
    throw precice::Error("You called the API function \"writeGradientData\", which is part of the experimental API. "
                         "You may unlock the full API by specifying <precice-configuration experimental=\"true\" ... > in the configuration. Please be aware that experimental features may change in any future version (even minor or bugfix).");
  }
  if (_impl->runtimeState.finalized) {
    throw precice::Error("writeGradientData() cannot be called after finalize().");
  }
  if (!_impl->runtimeState.initialized) {
    throw precice::Error("initialize() has to be called before writeGradientData().");
  }
  // Check that gradient data is actually enabled for this data
  if (!requiresGradientDataFor(meshName, dataName)) {
    std::string dataNameStr = toStr(dataName);
    throw precice::Error(precice::utils::format_or_error(
        "Data \"{}\" has no gradient values available. Please set the gradient flag to true under the data attribute in the configuration file.",
        dataNameStr));
  }
  if (ids.empty() && gradients.empty()) {
    return;
  }
  const auto meshDimensions = getMeshDimensions(meshName);
  const auto dataDimensions = getDataDimensions(meshName, dataName);
  const auto expectedArity  = static_cast<std::size_t>(meshDimensions * dataDimensions);
  if (gradients.size() != ids.size() * expectedArity) {
    throw precice::Error(precice::utils::format_or_error(
        "writeGradientData() was called with {} vertices and {} gradient values, but needs {} values ({} x {} x {}).",
        ids.size(), gradients.size(), ids.size() * expectedArity, ids.size(), meshDimensions, dataDimensions));
  }
  // no-op
}

// Profiling

void Participant::startProfilingSection(::precice::string_view sectionName)
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  std::string                           name = toStr(sectionName);
  if (name.find('/') != std::string::npos) {
    throw precice::Error(precice::utils::format_or_error(
        "The provided section name \"{}\" may not contain a forward-slash \"/\"", name));
  }
  ++_impl->runtimeState.activeProfilingSections;
}

void Participant::stopLastProfilingSection()
{
  std::lock_guard<std::recursive_mutex> lk(_impl->mtx);
  if (_impl->runtimeState.activeProfilingSections == 0) {
    throw precice::Error("There is no user-started event to stop.");
  }
  --_impl->runtimeState.activeProfilingSections;
}

// Missing Symbols

extern char const *const versionInformation __attribute__((visibility("default"))) = mockVersionInformation.c_str(); // Tooling API

namespace tooling {

void printConfigReference(std::ostream &out, ConfigReferenceType reftype)
{
  out << "preCICE mock config reference; type=" << static_cast<int>(reftype) << "\n";
}

void checkConfiguration(const std::string &filename, const std::string &participant, int size)
{
  throw precice::Error(precice::utils::format_or_error(
      "[precice-mock]  precice mock: checkConfiguration called for file '{}' participant '{}' size {}",
      filename, participant, size));
}

} // namespace tooling

std::string getVersionInformation()
{
  return mockVersionInformation;
}

} // namespace precice
