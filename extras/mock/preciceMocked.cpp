#include <chrono>
#include <fstream>
#include <future>
#include <iostream> //cpp11
#include <memory>   //make_unique from cpp14 rest cpp11
#include <mutex>    //cpp11
#include <random>
#include <string>
#include <thread>

#include <precice/Tooling.hpp>
#include <precice/Types.hpp>
#include <precice/span.hpp>
// project dependencies, used to remove reliance an std types maybe not available on older versions of cpp or outside libs

// --figure out why compilation was weird (weird issue when creating the current git structure merged old version of makefile incorrectly)
// --Ltrans warning is not easily avoidable (has to do with the way linking is handled in the entire project)
// --only happens on build from a clean state

// Add dummy data exchange
// user specifies timeseries of data from file
// run an open foam example against the mock (flap has this data alt. flow over heated plate)

// dont use random data instead use linear transform

// --add pre-commit

// ssh link for git

// API declaration
namespace precice {

// Using definition from Participant.hpp
using string_view = ::precice::span<const char>;

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

// Implementation
namespace precice {

namespace impl {

class FileDoubleReader {
public:
  explicit FileDoubleReader(const std::string &path = std::string()) : path_(path)
  {
    if (!path_.empty())
      open(path_);
  }

  bool open(const std::string &path)
  {
    std::lock_guard<std::mutex> lk(mtx_);
    path_ = path;
    ifs_.close();
    ifs_.clear();
    ifs_.open(path_, std::ios::in);
    return ifs_.is_open();
  }

  bool isOpen() const
  {
    return ifs_.is_open();
  }

  bool readNext(double &value)
  {
    std::lock_guard<std::mutex> lk(mtx_);
    if (!ifs_.is_open()) {
      if (path_.empty() || !open(path_))
        return false;
    }
    for (;;) {
      if (ifs_ >> value) {
        return true;
      }
      if (ifs_.eof()) {
        ifs_.clear();
        ifs_.seekg(0);
        if (ifs_.peek() == EOF)
          return false;
        continue;
      }
      return false;
    }
  }

private:
  std::string           path_;
  mutable std::ifstream ifs_;
  mutable std::mutex    mtx_;
};

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
    seed = static_cast<uint32_t>(rank) ^ 0x9e3779b9u;
  }

  ~ParticipantImpl() = default;

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
  std::string         promptedFilepath;
  bool                filepathSet = false;

  std::unique_ptr<FileDoubleReader> fileReader;
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
  if (!_impl)
    return;
  std::lock_guard<std::mutex> lk(_impl->mtx);

  if (!_impl->filepathSet) {
    std::cout << "Enter filepath for data series (10s timeout): " << std::flush;

    std::atomic<bool>         timedOut{false};
    std::promise<std::string> p;
    auto                      fut = p.get_future();

    std::thread t([&p, &timedOut]() {
      std::string input;
      if (std::getline(std::cin, input)) {
        if (!timedOut.load()) {
          p.set_value(input);
        }
      } else {
        if (!timedOut.load()) {
          p.set_value(std::string());
        }
      }
    });

    if (fut.wait_for(std::chrono::seconds(10)) == std::future_status::ready) {
      _impl->promptedFilepath = fut.get();
      _impl->useFile          = !_impl->promptedFilepath.empty();
      std::cout << "\nReceived filepath: " << _impl->promptedFilepath << std::endl;
      if (_impl->useFile) {
        _impl->fileReader = std::make_unique<impl::FileDoubleReader>(_impl->promptedFilepath);
        if (!_impl->fileReader->isOpen()) {
          std::cerr << "precice mock: failed to open file '" << _impl->promptedFilepath << "'. Falling back to random data.\n";
          _impl->useFile = false;
          _impl->fileReader.reset();
        }
      }
      if (t.joinable())
        t.join();
    } else {
      timedOut = true;
      std::cout << "\nNo input within 10 seconds. Proceeding without filepath.\n";
      if (t.joinable())
        t.detach();
    }
    _impl->filepathSet = true;
  }

  _impl->initialized     = true;
  _impl->couplingOngoing = true;
  _impl->maxTimeStep     = 1.0;
  _impl->currentStep     = 0;
  _impl->currentTime     = 0.0;
}

void Participant::advance(double computedTimeStepSize)
{
  if (!_impl)
    return;
  std::lock_guard<std::mutex> lk(_impl->mtx);
  (void) computedTimeStepSize;
  if (_impl->maxTimeStep == _impl->currentStep) {
    _impl->couplingOngoing = false;
    return;
  }

  _impl->currentStep += 1;
  _impl->currentTime += computedTimeStepSize;

  if (computedTimeStepSize > 0.0) {
    _impl->maxTimeStep = std::max(_impl->maxTimeStep, computedTimeStepSize);
  }
}

void Participant::finalize()
{
  if (!_impl)
    return;
  std::lock_guard<std::mutex> lk(_impl->mtx);
  _impl->initialized     = false;
  _impl->couplingOngoing = false;
}

// Implicit coupling

bool Participant::requiresWritingCheckpoint()
{
  return false;
}

bool Participant::requiresReadingCheckpoint()
{
  return false;
}

// Status queries

int Participant::getMeshDimensions(::precice::string_view /*meshName*/) const
{
  // mock in 3D
  return 3;
}

int Participant::getDataDimensions(::precice::string_view /*meshName*/,
                                   ::precice::string_view /*dataName*/) const
{
  // mock as scalar scalar
  return 1;
}

bool Participant::isCouplingOngoing() const
{
  if (!_impl)
    return false;
  return _impl->couplingOngoing;
}

bool Participant::isTimeWindowComplete() const
{
  return !isCouplingOngoing();
}

double Participant::getMaxTimeStepSize() const
{
  if (!_impl)
    return 0.0;
  return _impl->maxTimeStep;
}

// Mesh access

bool Participant::requiresMeshConnectivityFor(::precice::string_view /*meshName*/) const
{
  return false;
}

void Participant::resetMesh(::precice::string_view /*meshName*/)
{
  // no-op
}

VertexID Participant::setMeshVertex(
    ::precice::string_view /*meshName*/,
    ::precice::span<const double> /*position*/)
{
  return VertexID{};
}

int Participant::getMeshVertexSize(::precice::string_view /*meshName*/) const
{
  return 0;
}

void Participant::setMeshVertices(
    ::precice::string_view /*meshName*/,
    ::precice::span<const double> /*coordinates*/,
    ::precice::span<VertexID> /*ids*/)
{
  // no-op
}

void Participant::setMeshEdge(
    ::precice::string_view /*meshName*/,
    VertexID /*first*/,
    VertexID /*second*/)
{
  // no-op
}

void Participant::setMeshEdges(
    ::precice::string_view /*meshName*/,
    ::precice::span<const VertexID> /*ids*/)
{
  // no-op
}

void Participant::setMeshTriangle(
    ::precice::string_view /*meshName*/,
    VertexID /*first*/,
    VertexID /*second*/,
    VertexID /*third*/)
{
  // no-op
}

void Participant::setMeshTriangles(
    ::precice::string_view /*meshName*/,
    ::precice::span<const VertexID> /*ids*/)
{
  // no-op
}

void Participant::setMeshQuad(
    ::precice::string_view /*meshName*/,
    VertexID /*first*/,
    VertexID /*second*/,
    VertexID /*third*/,
    VertexID /*fourth*/)
{
  // no-op
}

void Participant::setMeshQuads(
    ::precice::string_view /*meshName*/,
    ::precice::span<const VertexID> /*ids*/)
{
  // no-op
}

void Participant::setMeshTetrahedron(
    ::precice::string_view /*meshName*/,
    VertexID /*first*/,
    VertexID /*second*/,
    VertexID /*third*/,
    VertexID /*fourth*/)
{
  // no-op
}

void Participant::setMeshTetrahedra(
    ::precice::string_view /*meshName*/,
    ::precice::span<const VertexID> /*ids*/)
{
  // no-op
}

// Data access

bool Participant::requiresInitialData()
{
  return false;
}

void Participant::writeData(
    ::precice::string_view /*meshName*/,
    ::precice::string_view /*dataName*/,
    ::precice::span<const VertexID> /*ids*/,
    ::precice::span<const double> /*values*/)
{
  // no-op
}

void Participant::readData(
    ::precice::string_view /*meshName*/,
    ::precice::string_view /*dataName*/,
    ::precice::span<const VertexID> /*ids*/,
    double /*relativeReadTime*/,
    ::precice::span<double> values) const
{
  if (!_impl)
    return;
  std::lock_guard<std::mutex> lk(_impl->mtx);
  const std::size_t           n = values.size();

  if (_impl->useFile && _impl->fileReader) // return data from file
  {
    for (std::size_t i = 0; i < n; ++i) {
      double val = 0.0;
      if (_impl->fileReader->readNext(val)) {
        values[i] = val;
      } else {
        values[i] = 0.0;
      }
    }
  } else // use random data
  {
    std::mt19937                           gen(static_cast<uint32_t>(_impl->seed + static_cast<uint32_t>(_impl->currentStep)));
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (std::size_t i = 0; i < n; ++i) {
      values[i] = dist(gen);
    }
  }
}

// Just-in-time mapping (experimental)

void Participant::writeAndMapData(
    ::precice::string_view /*meshName*/,
    ::precice::string_view /*dataName*/,
    ::precice::span<const double> /*coordinates*/,
    ::precice::span<const double> /*values*/)
{
  // no-op
}

void Participant::mapAndReadData(
    ::precice::string_view /*meshName*/,
    ::precice::string_view /*dataName*/,
    ::precice::span<const double> /*coordinates*/,
    double /*relativeReadTime*/,
    ::precice::span<double> /*values*/) const
{
  // no-op
}

// Direct access

void Participant::setMeshAccessRegion(
    ::precice::string_view /*meshName*/,
    ::precice::span<const double> /*boundingBox*/) const
{
  // no-op
}

void Participant::getMeshVertexIDsAndCoordinates(
    ::precice::string_view /*meshName*/,
    ::precice::span<VertexID> /*ids*/,
    ::precice::span<double> /*coordinates*/) const
{
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
    ::precice::span<const VertexID> /*ids*/,
    ::precice::span<const double> /*gradients*/)
{
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
