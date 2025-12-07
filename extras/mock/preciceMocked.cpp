#include <chrono>
#include <fstream>
#include <future>
#include <iostream>
#include <memory>
#include <mutex>
#include <random>
#include <regex>
#include <sstream>
#include <string>
#include <thread>

#include <precice/Tooling.hpp>
#include <precice/Types.hpp>
#include <precice/span.hpp>
#include "precice/Exceptions.hpp"
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
  VertexID setMeshVertex(::precice::string_view        meshName,
                         ::precice::span<const double> position);
  int      getMeshVertexSize(::precice::string_view meshName) const;
  void     setMeshVertices(::precice::string_view        meshName,
                           ::precice::span<const double> coordinates,
                           ::precice::span<VertexID>     ids);
  void     setMeshEdge(::precice::string_view meshName, VertexID first, VertexID second);
  void     setMeshEdges(::precice::string_view          meshName,
                        ::precice::span<const VertexID> ids);
  void     setMeshTriangle(::precice::string_view meshName, VertexID first, VertexID second, VertexID third);
  void     setMeshTriangles(::precice::string_view          meshName,
                            ::precice::span<const VertexID> ids);
  void     setMeshQuad(::precice::string_view meshName, VertexID first, VertexID second, VertexID third, VertexID fourth);
  void     setMeshQuads(::precice::string_view          meshName,
                        ::precice::span<const VertexID> ids);
  void     setMeshTetrahedron(::precice::string_view meshName, VertexID first, VertexID second, VertexID third, VertexID fourth);
  void     setMeshTetrahedra(::precice::string_view          meshName,
                             ::precice::span<const VertexID> ids);

  // Data access
  bool requiresInitialData();
  void writeData(::precice::string_view          meshName,
                 ::precice::string_view          dataName,
                 ::precice::span<const VertexID> ids,
                 ::precice::span<const double>   values);
  void readData(::precice::string_view          meshName,
                ::precice::string_view          dataName,
                ::precice::span<const VertexID> ids,
                double                          relativeReadTime,
                ::precice::span<double>         values) const;

  // Just-in-time mapping (experimental)
  void writeAndMapData(::precice::string_view        meshName,
                       ::precice::string_view        dataName,
                       ::precice::span<const double> coordinates,
                       ::precice::span<const double> values);
  void mapAndReadData(::precice::string_view        meshName,
                      ::precice::string_view        dataName,
                      ::precice::span<const double> coordinates,
                      double                        relativeReadTime,
                      ::precice::span<double>       values) const;

  // Direct access
  void setMeshAccessRegion(::precice::string_view        meshName,
                           ::precice::span<const double> boundingBox) const;
  void getMeshVertexIDsAndCoordinates(::precice::string_view    meshName,
                                      ::precice::span<VertexID> ids,
                                      ::precice::span<double>   coordinates) const;

  // Gradient data (experimental)
  bool requiresGradientDataFor(::precice::string_view meshName,
                               ::precice::string_view dataName) const;
  void writeGradientData(::precice::string_view          meshName,
                         ::precice::string_view          dataName,
                         ::precice::span<const VertexID> ids,
                         ::precice::span<const double>   gradients);

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
  Currently only basic checks using assertions are performed.
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
                  void                  *communicator = nullptr)
      : name(std::string(participantName.data(), participantName.size())),
        config(std::string(configurationFileName.data(), configurationFileName.size())),
        rank(solverProcessIndex),
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
    int         dimensions = 1; // scalar by default
    bool        isRead     = false;
    bool        isWrite    = false;
  };

  struct ConfigData {
    std::map<std::string, MeshInfo> meshes;
    std::vector<DataInfo>           dataItems;
    bool                            isImplicitCoupling = false;
    int                             maxIterations      = 1;
    int                             currentIteration   = 0;
    bool                            iterationConverged = false;
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
    DataMode                              defaultMode = DataMode::Buffer;
  };

  // Minimal internal state
  std::string name;
  std::string config;
  int         rank            = 0;
  int         size            = 1;
  void       *comm            = nullptr;
  bool        initialized     = false;
  bool        finalized       = false;
  bool        couplingOngoing = false;
  double      maxTimeStep     = 0.0;
  std::mutex  mtx;

  std::vector<double> dataSeries;
  bool                useFile     = false;
  uint32_t            seed        = 0;
  std::size_t         currentStep = 0;
  double              currentTime = 0.0;

  ConfigData configData;
  bool       configParsed = false;

  MockConfig mockConfig;
  bool       mockConfigParsed = false;

  // Write data buffers: Key is "meshName:dataName"
  std::map<std::string, std::vector<double>> writeBuffers;

  // std::unique_ptr<FileDoubleReader> fileReader;

  void parseConfig()
  {
    if (configParsed)
      return;

    try {
      std::ifstream file(config);
      if (!file.is_open()) {
        throw precice::Error(precice::utils::format_or_error(
            "Could not open configuration file '{}'", config));
      }

      std::string line;
      std::string currentMesh, currentData;
      bool        inParticipant    = false;
      bool        foundParticipant = false;
      int         meshDims         = 3;

      // First pass: collect data type information (scalar vs vector)
      std::map<std::string, int> dataTypeDimensions;
      std::ifstream              fileFirst(config);
      while (std::getline(fileFirst, line)) {
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);

        if (line.find("<data:scalar") != std::string::npos && line.find("name=\"") != std::string::npos) {
          size_t      namePos          = line.find("name=\"") + 6;
          size_t      nameEnd          = line.find("\"", namePos);
          std::string dataName         = line.substr(namePos, nameEnd - namePos);
          dataTypeDimensions[dataName] = 1;
        } else if (line.find("<data:vector") != std::string::npos && line.find("name=\"") != std::string::npos) {
          size_t      namePos          = line.find("name=\"") + 6;
          size_t      nameEnd          = line.find("\"", namePos);
          std::string dataName         = line.substr(namePos, nameEnd - namePos);
          dataTypeDimensions[dataName] = 3; // Vectors are 3D
        }
      }
      fileFirst.close();

      // Second pass: parse participant section
      file.clear();
      file.seekg(0);
      while (std::getline(file, line)) {
        // Trim whitespace
        line.erase(0, line.find_first_not_of(" \t\r\n"));
        line.erase(line.find_last_not_of(" \t\r\n") + 1);

        // Check for participant tag
        if (line.find("<participant") != std::string::npos && line.find("name=\"") != std::string::npos) {
          size_t      namePos         = line.find("name=\"") + 6;
          size_t      nameEnd         = line.find("\"", namePos);
          std::string participantName = line.substr(namePos, nameEnd - namePos);

          if (participantName == name) {
            inParticipant    = true;
            foundParticipant = true;
          }
        } else if (line.find("</participant>") != std::string::npos) {
          inParticipant = false;
        }

        if (!inParticipant)
          continue;

        // Parse mesh tags
        if (line.find("<mesh") != std::string::npos && line.find("name=\"") != std::string::npos) {
          size_t namePos = line.find("name=\"") + 6;
          size_t nameEnd = line.find("\"", namePos);
          currentMesh    = line.substr(namePos, nameEnd - namePos);

          size_t dimsPos = line.find("dimensions=\"");
          if (dimsPos != std::string::npos) {
            dimsPos += 12;
            size_t dimsEnd = line.find("\"", dimsPos);
            meshDims       = std::stoi(line.substr(dimsPos, dimsEnd - dimsPos));
          }

          MeshInfo meshInfo;
          meshInfo.name                  = currentMesh;
          meshInfo.dimensions            = meshDims;
          configData.meshes[currentMesh] = meshInfo;
        }

        // Parse provide-mesh
        if (line.find("<provide-mesh") != std::string::npos && line.find("name=\"") != std::string::npos) {
          size_t      namePos  = line.find("name=\"") + 6;
          size_t      nameEnd  = line.find("\"", namePos);
          std::string meshName = line.substr(namePos, nameEnd - namePos);

          auto it = configData.meshes.find(meshName);
          if (it != configData.meshes.end()) {
            it->second.provided = true;
          }
        }

        // Parse receive-mesh
        if (line.find("<receive-mesh") != std::string::npos && line.find("name=\"") != std::string::npos) {
          size_t      namePos  = line.find("name=\"") + 6;
          size_t      nameEnd  = line.find("\"", namePos);
          std::string meshName = line.substr(namePos, nameEnd - namePos);

          auto it = configData.meshes.find(meshName);
          if (it != configData.meshes.end()) {
            it->second.received = true;
          }
        }

        // Parse write-data
        if (line.find("<write-data") != std::string::npos) {
          size_t namePos = line.find("name=\"");
          size_t meshPos = line.find("mesh=\"");

          if (namePos != std::string::npos && meshPos != std::string::npos) {
            namePos += 6;
            size_t      nameEnd  = line.find("\"", namePos);
            std::string dataName = line.substr(namePos, nameEnd - namePos);

            meshPos += 6;
            size_t      meshEnd  = line.find("\"", meshPos);
            std::string meshName = line.substr(meshPos, meshEnd - meshPos);

            DataInfo dataInfo;
            dataInfo.name     = dataName;
            dataInfo.meshName = meshName;
            dataInfo.isWrite  = true;
            // Look up dimensionality from data type declarations
            auto dimIt          = dataTypeDimensions.find(dataName);
            dataInfo.dimensions = (dimIt != dataTypeDimensions.end()) ? dimIt->second : 1;
            configData.dataItems.push_back(dataInfo);
          }
        }

        // Parse read-data
        if (line.find("<read-data") != std::string::npos) {
          size_t namePos = line.find("name=\"");
          size_t meshPos = line.find("mesh=\"");

          if (namePos != std::string::npos && meshPos != std::string::npos) {
            namePos += 6;
            size_t      nameEnd  = line.find("\"", namePos);
            std::string dataName = line.substr(namePos, nameEnd - namePos);

            meshPos += 6;
            size_t      meshEnd  = line.find("\"", meshPos);
            std::string meshName = line.substr(meshPos, meshEnd - meshPos);

            DataInfo dataInfo;
            dataInfo.name     = dataName;
            dataInfo.meshName = meshName;
            dataInfo.isRead   = true;
            // Look up dimensionality from data type declarations
            auto dimIt          = dataTypeDimensions.find(dataName);
            dataInfo.dimensions = (dimIt != dataTypeDimensions.end()) ? dimIt->second : 1;
            configData.dataItems.push_back(dataInfo);
          }
        }

        // Detect implicit coupling
        if (line.find("<coupling-scheme:serial-implicit") != std::string::npos) {
          configData.isImplicitCoupling = true;
        }
        if (line.find("<coupling-scheme:parallel-implicit") != std::string::npos) {
          configData.isImplicitCoupling = true;
        }
      }

      if (!foundParticipant) {
        throw precice::Error(precice::utils::format_or_error(
            "Participant '{}' not found in configuration '{}'", name, config));
      }

      // Set defaults for implicit coupling
      if (configData.isImplicitCoupling) {
        configData.maxIterations = 5;
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

  void parseMockConfig()
  {
    if (mockConfigParsed)
      return;

    try {
      // Derive mock-config path from precice config path
      std::string mockConfigPath = config;
      size_t      lastSlash      = mockConfigPath.find_last_of("/\\");
      size_t      lastDot        = mockConfigPath.find_last_of('.');

      if (lastSlash != std::string::npos) {
        // Has directory path
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
        return;
      }
      mockFile.close();

      // Parse mock config using simple XML parser
      // Format: <mock-config>
      //           <data-mapping mesh="MeshName" data="DataName" mode="buffer|random|scaled">
      //             <scalar-multiplier value="2.0" /> or
      //             <vector-multiplier values="1.0;2.0;3.0" />
      //           </data-mapping>
      //         </mock-config>

      std::ifstream       file(mockConfigPath);
      std::string         line;
      std::string         currentMesh, currentData;
      DataMode            currentMode   = DataMode::Buffer;
      double              currentScalar = 1.0;
      std::vector<double> currentVector;

      while (std::getline(file, line)) {
        // Simple parsing - look for key attributes
        if (line.find("<data-mapping") != std::string::npos) {
          // Extract mesh
          size_t meshPos = line.find("mesh=\"");
          if (meshPos != std::string::npos) {
            meshPos += 6;
            size_t meshEnd = line.find("\"", meshPos);
            currentMesh    = line.substr(meshPos, meshEnd - meshPos);
          }

          // Extract data
          size_t dataPos = line.find("data=\"");
          if (dataPos != std::string::npos) {
            dataPos += 6;
            size_t dataEnd = line.find("\"", dataPos);
            currentData    = line.substr(dataPos, dataEnd - dataPos);
          }

          // Extract mode
          size_t modePos = line.find("mode=\"");
          if (modePos != std::string::npos) {
            modePos += 6;
            size_t      modeEnd = line.find("\"", modePos);
            std::string modeStr = line.substr(modePos, modeEnd - modePos);
            if (modeStr == "random") {
              currentMode = DataMode::Random;
            } else if (modeStr == "scaled") {
              currentMode = DataMode::ScaledBuffer;
            } else {
              currentMode = DataMode::Buffer;
            }
          }

          currentScalar = 1.0;
          currentVector.clear();
        } else if (line.find("<scalar-multiplier") != std::string::npos) {
          size_t valPos = line.find("value=\"");
          if (valPos != std::string::npos) {
            valPos += 7;
            size_t valEnd = line.find("\"", valPos);
            currentScalar = std::stod(line.substr(valPos, valEnd - valPos));
          }
        } else if (line.find("<vector-multiplier") != std::string::npos) {
          size_t valPos = line.find("values=\"");
          if (valPos != std::string::npos) {
            valPos += 8;
            size_t      valEnd    = line.find("\"", valPos);
            std::string valuesStr = line.substr(valPos, valEnd - valPos);

            // Parse semicolon-separated values
            size_t pos = 0;
            while (pos < valuesStr.length()) {
              size_t nextSemi = valuesStr.find(';', pos);
              if (nextSemi == std::string::npos)
                nextSemi = valuesStr.length();
              std::string valStr = valuesStr.substr(pos, nextSemi - pos);
              currentVector.push_back(std::stod(valStr));
              pos = nextSemi + 1;
            }
          }
        } else if (line.find("</data-mapping>") != std::string::npos) {
          if (!currentMesh.empty() && !currentData.empty()) {
            std::string    key = currentMesh + ":" + currentData;
            MockDataConfig config;
            config.dataName             = currentData;
            config.meshName             = currentMesh;
            config.mode                 = currentMode;
            config.scalarMultiplier     = currentScalar;
            config.vectorMultiplier     = currentVector;
            mockConfig.dataConfigs[key] = config;
          }
        }
      }

      mockConfigParsed = true;
    } catch (const std::exception &e) {
      std::cerr << "Warning: Could not parse mock configuration file: " << e.what()
                << ". Using default buffer mode." << std::endl;
      mockConfig.defaultMode = DataMode::Buffer;
      mockConfigParsed       = true;
    }
  }
};

} // namespace impl

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
}

Participant::~Participant() = default;

// Steering methods

void Participant::initialize()
{
  if (!_impl) {
    throw precice::Error("Participant implementation missing in initialize().");
  }
  std::lock_guard<std::mutex> lk(_impl->mtx);
  if (_impl->finalized) {
    throw precice::Error("initialize() cannot be called after finalize().");
  }
  if (_impl->initialized) {
    throw precice::Error("initialize() may only be called once.");
  }

  // Parse configuration
  _impl->parseConfig();

  // Parse mock configuration
  _impl->parseMockConfig();

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
  if (!_impl) {
    throw precice::Error("Participant implementation missing in advance().");
  }
  std::lock_guard<std::mutex> lk(_impl->mtx);
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
  if (!_impl) {
    throw precice::Error("Participant implementation missing in finalize().");
  }
  std::lock_guard<std::mutex> lk(_impl->mtx);
  if (_impl->finalized) {
    throw precice::Error("finalize() may only be called once.");
  }
  _impl->initialized     = false;
  _impl->couplingOngoing = false;
}

// Implicit coupling

bool Participant::requiresWritingCheckpoint()
{
  if (!_impl) {
    throw precice::Error("Participant implementation missing in requiresWritingCheckpoint().");
  }
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
  if (!_impl) {
    throw precice::Error("Participant implementation missing in requiresReadingCheckpoint().");
  }
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
  if (!_impl) {
    throw precice::Error("Participant implementation missing in getMeshDimensions().");
  }

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

  // Fallback: mock in 3D
  return 3;
}

int Participant::getDataDimensions(::precice::string_view meshName,
                                   ::precice::string_view dataName) const
{
  if (!_impl) {
    throw precice::Error("Participant implementation missing in getDataDimensions().");
  }

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
  if (!_impl) {
    throw precice::Error("Participant implementation missing in isCouplingOngoing().");
  }
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
  if (!_impl) {
    throw precice::Error("Participant implementation missing in isTimeWindowComplete().");
  }
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
  if (!_impl) {
    throw precice::Error("Participant implementation missing in getMaxTimeStepSize().");
  }
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

void Participant::resetMesh(::precice::string_view /*meshName*/)
{
  if (!_impl) {
    throw precice::Error("Participant implementation missing in resetMesh().");
  }
  if (!_impl->initialized) {
    throw precice::Error("initialize() has to be called before resetMesh().");
  }
  // no-op
}

VertexID Participant::setMeshVertex(
    ::precice::string_view /*meshName*/,
    ::precice::span<const double> position)
{
  if (!_impl) {
    throw precice::Error("Participant implementation missing in setMeshVertex().");
  }
  if (_impl->initialized) {
    throw precice::Error("setMeshVertex() cannot be called after initialize(). Mesh modification is only allowed before calling initialize().");
  }
  // Expect a 3D position for a single vertex in this mock
  if (position.size() != 3) {
    throw precice::Error(precice::utils::format_or_error("setMeshVertex() was called with {} coordinates, but expects 3 coordinates in 3D.", position.size()));
  }
  return VertexID{};
}

int Participant::getMeshVertexSize(::precice::string_view /*meshName*/) const
{
  return 0;
}

void Participant::setMeshVertices(
    ::precice::string_view /*meshName*/,
    ::precice::span<const double> coordinates,
    ::precice::span<VertexID>     ids)
{
  if (!_impl) {
    throw precice::Error("Participant implementation missing in setMeshVertices().");
  }
  if (_impl->initialized) {
    throw precice::Error("setMeshVertices() cannot be called after initialize(). Mesh modification is only allowed before calling initialize().");
  }
  // Expect coordinates as 3 * N entries for N VertexIDs
  if (coordinates.size() != ids.size() * 3) {
    throw precice::Error(precice::utils::format_or_error("setMeshVertices() was called with {} vertices and {} coordinates ({}D), but needs {} coordinates ({} x 3).", ids.size(), coordinates.size(), 3, ids.size() * 3, ids.size()));
  }
  // no-op
}

void Participant::setMeshEdge(
    ::precice::string_view /*meshName*/,
    VertexID first,
    VertexID second)
{
  if (!_impl) {
    throw precice::Error("Participant implementation missing in setMeshEdge().");
  }
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
    ::precice::span<const VertexID> ids)
{
  if (!_impl) {
    throw precice::Error("Participant implementation missing in setMeshEdges().");
  }
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
  if (!_impl) {
    throw precice::Error("Participant implementation missing in setMeshTriangle().");
  }
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
    ::precice::span<const VertexID> ids)
{
  if (!_impl) {
    throw precice::Error("Participant implementation missing in setMeshTriangles().");
  }
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
  if (!_impl) {
    throw precice::Error("Participant implementation missing in setMeshQuad().");
  }
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
    ::precice::span<const VertexID> ids)
{
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
  if (!_impl) {
    throw precice::Error("Participant implementation missing in setMeshTetrahedron().");
  }
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
    ::precice::span<const VertexID> ids)
{
  if (!_impl) {
    throw precice::Error("Participant implementation missing in setMeshTetrahedra().");
  }
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
  if (!_impl) {
    throw precice::Error("Participant implementation missing in requiresInitialData().");
  }
  if (_impl->initialized) {
    throw precice::Error("requiresInitialData() has to be called before initialize().");
  }
  return false;
}

void Participant::writeData(
    ::precice::string_view          meshName,
    ::precice::string_view          dataName,
    ::precice::span<const VertexID> ids,
    ::precice::span<const double>   values)
{
  if (!_impl) {
    throw precice::Error("Participant implementation missing in writeData().");
  }
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
  std::string                 key = meshNameStr + ":" + dataNameStr;
  std::lock_guard<std::mutex> lk(_impl->mtx);
  _impl->writeBuffers[key].assign(values.begin(), values.end());
}

void Participant::readData(
    ::precice::string_view          meshName,
    ::precice::string_view          dataName,
    ::precice::span<const VertexID> ids,
    double /*relativeReadTime*/,
    ::precice::span<double> values) const
{
  if (!_impl) {
    throw precice::Error("Participant implementation missing in readData().");
  }
  std::lock_guard<std::mutex> lk(_impl->mtx);
  const std::size_t           n = values.size();
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

  auto mockIt = _impl->mockConfig.dataConfigs.find(key);
  if (mockIt != _impl->mockConfig.dataConfigs.end()) {
    mode           = mockIt->second.mode;
    mockDataConfig = &mockIt->second;
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
      // Copy as much as we can from the buffer
      size_t copySize = std::min(n, buffer.size());
      for (std::size_t i = 0; i < copySize; ++i) {
        values[i] = buffer[i];
      }
      // Fill remaining with zeros if buffer is smaller
      for (std::size_t i = copySize; i < n; ++i) {
        values[i] = 0.0;
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
      const auto &buffer   = bufferIt->second;
      size_t      copySize = std::min(n, buffer.size());

      if (mockDataConfig && !mockDataConfig->vectorMultiplier.empty()) {
        // Element-wise multiplication with vector
        const auto &multiplier = mockDataConfig->vectorMultiplier;
        for (std::size_t i = 0; i < copySize; ++i) {
          values[i] = buffer[i] * multiplier[i % multiplier.size()];
        }
      } else {
        // Scalar multiplication
        double scalar = mockDataConfig ? mockDataConfig->scalarMultiplier : 1.0;
        for (std::size_t i = 0; i < copySize; ++i) {
          values[i] = buffer[i] * scalar;
        }
      }

      // Fill remaining with zeros
      for (std::size_t i = copySize; i < n; ++i) {
        values[i] = 0.0;
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
    ::precice::span<const double> coordinates,
    ::precice::span<const double> values)
{
  if (!_impl) {
    throw precice::Error("Participant implementation missing in writeAndMapData().");
  }
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
    ::precice::span<const double> coordinates,
    double /*relativeReadTime*/,
    ::precice::span<double> values) const
{
  if (!_impl) {
    throw precice::Error("Participant implementation missing in mapAndReadData().");
  }
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
    ::precice::span<const double> boundingBox) const
{
  if (!_impl) {
    throw precice::Error("Participant implementation missing in setMeshAccessRegion().");
  }
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
    ::precice::span<VertexID> ids,
    ::precice::span<double>   coordinates) const
{
  if (!_impl) {
    throw precice::Error("Participant implementation missing in getMeshVertexIDsAndCoordinates().");
  }
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
    ::precice::span<const VertexID> ids,
    ::precice::span<const double>   gradients)
{
  if (!_impl) {
    throw precice::Error("Participant implementation missing in writeGradientData().");
  }
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
  if (!_impl) {
    throw precice::Error("Participant implementation missing in startProfilingSection().");
  }
  // no-op
}

void Participant::stopLastProfilingSection()
{
  if (!_impl) {
    throw precice::Error("Participant implementation missing in stopLastProfilingSection().");
  }
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
