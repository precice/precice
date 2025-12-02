#include <cassert>
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

/*class FileDoubleReader {
public:
  explicit FileDoubleReader(const std::string &path = std::string()) : path_(path)
  {
    if (!path_.empty()) //rewrite using assertion
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
};*/

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
    assert(solverProcessSize > 0 && "solverProcessSize must be > 0");
    assert(solverProcessIndex >= 0 && solverProcessIndex < solverProcessSize && "solverProcessIndex must be in [0, solverProcessSize)");
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

  // std::unique_ptr<FileDoubleReader> fileReader;
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
  assert(_impl && "Participant implementation missing in initialize().");
  if (!_impl)
    return;
  std::lock_guard<std::mutex> lk(_impl->mtx);
  assert(!_impl->initialized && "initialize() called while participant is already initialized.");
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
      /*if (_impl->useFile) {
        _impl->fileReader = std::make_unique<impl::FileDoubleReader>(_impl->promptedFilepath);
        if (!_impl->fileReader->isOpen()) {
          std::cerr << "precice mock: failed to open file '" << _impl->promptedFilepath << "'. Falling back to random data.\n";
          _impl->useFile = false;
          _impl->fileReader.reset();
        }
      }*/
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
  assert(_impl && "Participant implementation missing in advance().");
  if (!_impl)
    return;
  std::lock_guard<std::mutex> lk(_impl->mtx);
  (void) computedTimeStepSize;
  assert(_impl->initialized && "advance() called before initialize().");
  assert(computedTimeStepSize >= 0.0 && "computedTimeStepSize must be non-negative.");
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
  assert(_impl && "Participant implementation missing in finalize().");
  if (!_impl)
    return;
  std::lock_guard<std::mutex> lk(_impl->mtx);
  assert(_impl->initialized && "finalize() called before initialize().");
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
  assert(_impl && "Participant implementation missing in isCouplingOngoing().");
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
  assert(_impl && "Participant implementation missing in getMaxTimeStepSize().");
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
    ::precice::span<const double> position)
{
  // Expect a 3D position for a single vertex in this mock
  assert(position.size() == 3 && "setMeshVertex: position must be of size 3 (3D) in the mock implementation.");
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
  // Expect coordinates as 3 * N entries for N VertexIDs
  assert(coordinates.size() == (ids.size() * 3) && "setMeshVertices: coordinates should have 3*ids.size() entries in the mock implementation.");
  // no-op
}

void Participant::setMeshEdge(
    ::precice::string_view /*meshName*/,
    VertexID first,
    VertexID second)
{
  // basic check: edge endpoints should not be identical
  assert(first != second && "setMeshEdge: first and second vertex ids should be different.");
  // no-op
}

void Participant::setMeshEdges(
    ::precice::string_view /*meshName*/,
    ::precice::span<const VertexID> ids)
{
  // expect pairs of vertex ids
  assert((ids.size() % 2 == 0) && "setMeshEdges: ids size must be divisible by 2.");
  // no-op
}

void Participant::setMeshTriangle(
    ::precice::string_view /*meshName*/,
    VertexID first,
    VertexID second,
    VertexID third)
{
  // triangle vertices should be distinct
  assert(first != second && first != third && second != third && "setMeshTriangle: triangle vertices must be distinct.");
  // no-op
}

void Participant::setMeshTriangles(
    ::precice::string_view /*meshName*/,
    ::precice::span<const VertexID> ids)
{
  // expect triples of ids
  assert((ids.size() % 3 == 0) && "setMeshTriangles: ids size must be divisible by 3.");
  // no-op
}

void Participant::setMeshQuad(
    ::precice::string_view /*meshName*/,
    VertexID first,
    VertexID second,
    VertexID third,
    VertexID fourth)
{
  // quad vertices should be distinct
  assert(first != second && first != third && first != fourth && second != third && second != fourth && third != fourth && "setMeshQuad: quad vertices should be distinct.");
  // no-op
}

void Participant::setMeshQuads(
    ::precice::string_view /*meshName*/,
    ::precice::span<const VertexID> ids)
{
  // expect groups of 4 vertex ids
  assert((ids.size() % 4 == 0) && "setMeshQuads: ids size must be divisible by 4.");
  // no-op
}

void Participant::setMeshTetrahedron(
    ::precice::string_view /*meshName*/,
    VertexID first,
    VertexID second,
    VertexID third,
    VertexID fourth)
{
  // tetrahedron vertices should be distinct
  assert(first != second && first != third && first != fourth && second != third && second != fourth && third != fourth && "setMeshTetrahedron: vertices must be distinct.");
  // no-op
}

void Participant::setMeshTetrahedra(
    ::precice::string_view /*meshName*/,
    ::precice::span<const VertexID> ids)
{
  // expect groups of 4 vertex ids
  assert((ids.size() % 4 == 0) && "setMeshTetrahedra: ids size must be divisible by 4.");
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
    ::precice::span<const VertexID> ids,
    ::precice::span<const double>   values)
{
  assert(_impl && "Participant implementation missing in writeData().");
  assert(_impl->initialized && "writeData() called before initialize().");
  // expect one value per id (scalar data in the mock)
  assert(values.size() == ids.size() && "writeData: values size must equal number of ids.");
  // no-op
}

void Participant::readData(
    ::precice::string_view /*meshName*/,
    ::precice::string_view /*dataName*/,
    ::precice::span<const VertexID> ids,
    double /*relativeReadTime*/,
    ::precice::span<double> values) const
{
  assert(_impl && "Participant implementation missing in readData().");
  if (!_impl)
    return;
  std::lock_guard<std::mutex> lk(_impl->mtx);
  const std::size_t           n = values.size();
  assert(_impl->initialized && "readData() called before initialize().");
  assert(n == ids.size() && "readData: values size must equal number of ids.");

  if (_impl->useFile /*&& _impl->fileReader*/) // return data from file
  {
    /*for (std::size_t i = 0; i < n; ++i) {
      double val = 0.0;
      if (_impl->fileReader->readNext(val)) {
        values[i] = val;
      } else {
        values[i] = 0.0;
      }
    }*/
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
    ::precice::span<const double> coordinates,
    ::precice::span<const double> values)
{
  assert(_impl && "Participant implementation missing in writeAndMapData().");
  assert(_impl->initialized && "writeAndMapData() called before initialize().");
  // coordinates are 3*N, values are N for scalar data
  assert(coordinates.size() % 3 == 0 && "writeAndMapData: coordinates size must be a multiple of 3.");
  assert(values.size() == coordinates.size() / 3 && "writeAndMapData: values size must equal number of coordinate points.");
  // no-op
}

void Participant::mapAndReadData(
    ::precice::string_view /*meshName*/,
    ::precice::string_view /*dataName*/,
    ::precice::span<const double> coordinates,
    double /*relativeReadTime*/,
    ::precice::span<double> values) const
{
  assert(_impl && "Participant implementation missing in mapAndReadData().");
  assert(_impl->initialized && "mapAndReadData() called before initialize().");
  assert(coordinates.size() % 3 == 0 && "mapAndReadData: coordinates size must be a multiple of 3.");
  assert(values.size() == coordinates.size() / 3 && "mapAndReadData: values size must equal number of coordinate points.");
  // no-op
}

// Direct access

void Participant::setMeshAccessRegion(
    ::precice::string_view /*meshName*/,
    ::precice::span<const double> boundingBox) const
{
  assert(_impl && "Participant implementation missing in setMeshAccessRegion().");
  assert(_impl->initialized && "setMeshAccessRegion() called before initialize().");
  // bounding box: min/max for each dimension -> 3D: 6 entries
  assert(boundingBox.size() == 6 && "setMeshAccessRegion: expected 6 entries for 3D bounding box (min/max per dimension).");
  // no-op
}

void Participant::getMeshVertexIDsAndCoordinates(
    ::precice::string_view /*meshName*/,
    ::precice::span<VertexID> ids,
    ::precice::span<double>   coordinates) const
{
  assert(_impl && "Participant implementation missing in getMeshVertexIDsAndCoordinates().");
  assert(_impl->initialized && "getMeshVertexIDsAndCoordinates() called before initialize().");
  // coordinates are 3 * ids.size()
  assert(coordinates.size() == ids.size() * 3 && "getMeshVertexIDsAndCoordinates: coordinates size must be 3 * ids.size().");
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
  assert(_impl && "Participant implementation missing in writeGradientData().");
  assert(_impl->initialized && "writeGradientData() called before initialize().");
  // gradients expected as 3 components per id for vector gradient in 3D
  assert(gradients.size() == ids.size() * 3 && "writeGradientData: gradients size must be 3 * ids.size() in this mock.");
  // no-op
}

// Profiling

void Participant::startProfilingSection(::precice::string_view /*sectionName*/)
{
  assert(_impl && "Participant implementation missing in startProfilingSection().");
  assert(_impl->initialized && "startProfilingSection() called before initialize().");
  // no-op
}

void Participant::stopLastProfilingSection()
{
  assert(_impl && "Participant implementation missing in stopLastProfilingSection().");
  assert(_impl->initialized && "stopLastProfilingSection() called before initialize().");
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
