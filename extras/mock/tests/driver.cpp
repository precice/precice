// Differential test driver: exercises the preCICE Participant API and logs
// observable behavior (return values / exception messages) with a "T> " prefix.
//
// Run the same binary against libprecice.so (as a coupled two-participant pair)
// and against libpreciceMocked.so (LD_PRELOAD, standalone) and diff the T> lines.
// See run.sh for the full scenario matrix.
#include <precice/precice.hpp>

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using precice::Participant;
using precice::VertexID;

static void logline(const std::string &s)
{
  std::cout << "T> " << s << std::endl;
}

#define LOG(expr)                 \
  do {                            \
    std::ostringstream _os;       \
    _os << expr;                  \
    logline(_os.str());           \
  } while (0)

static std::string fmtv(const std::vector<double> &v)
{
  std::ostringstream os;
  os.precision(6);
  os << std::fixed;
  for (size_t i = 0; i < v.size(); ++i) {
    if (i)
      os << ",";
    os << v[i];
  }
  return os.str();
}

template <typename F>
static void expectError(const std::string &what, F f)
{
  try {
    f();
    LOG("ERR " << what << ": NO-EXCEPTION");
  } catch (const std::exception &e) {
    LOG("ERR " << what << ": " << e.what());
  }
}

struct Setup {
  std::string me, other, myMesh, otherMesh, writeData, readData;
  int         writeDims, readDims;
};

static Setup makeSetup(const std::string &participant)
{
  // DataOne: scalar, written by SolverOne. DataTwo: vector (2D), written by SolverTwo.
  if (participant == "SolverOne") {
    return {"SolverOne", "SolverTwo", "MeshOne", "MeshTwo", "DataOne", "DataTwo", 1, 2};
  }
  return {"SolverTwo", "SolverOne", "MeshTwo", "MeshOne", "DataTwo", "DataOne", 2, 1};
}

static std::vector<VertexID> setupMesh(Participant &p, const Setup &s)
{
  // 3 vertices on a line, identical for both participants (identity NN mapping)
  std::vector<double>   coords = {0.0, 0.0, 1.0, 0.0, 2.0, 0.0};
  std::vector<VertexID> ids(3);
  p.setMeshVertices(s.myMesh, coords, ids);
  return ids;
}

static std::vector<double> makeValues(const Setup &s, int step)
{
  std::vector<double> v(3 * s.writeDims);
  for (size_t j = 0; j < v.size(); ++j) {
    v[j] = (step + 1) + 0.01 * static_cast<double>(j);
  }
  return v;
}

// Common query block after construction
static void logStaticQueries(Participant &p, const Setup &s)
{
  LOG("getMeshDimensions(" << s.myMesh << ")=" << p.getMeshDimensions(s.myMesh));
  LOG("getDataDimensions(write)=" << p.getDataDimensions(s.myMesh, s.writeData));
  LOG("getDataDimensions(read)=" << p.getDataDimensions(s.myMesh, s.readData));
  LOG("requiresMeshConnectivityFor=" << p.requiresMeshConnectivityFor(s.myMesh));
  LOG("requiresGradientDataFor(write)=" << p.requiresGradientDataFor(s.myMesh, s.writeData));
}

// ---------------------------------------------------------------------------
// Scenario: full lifecycle (works for explicit and implicit configs)
// fixedDt > 0: solver-prescribed step (for method="first-participant" configs)
static int runLifecycle(const std::string &name, const std::string &config, bool subcycle,
                        double fixedDt = -1.0)
{
  Setup       s = makeSetup(name);
  Participant p(name, config, 0, 1);

  logStaticQueries(p, s);

  auto ids = setupMesh(p, s);
  LOG("vertexIDs=" << ids[0] << "," << ids[1] << "," << ids[2]);
  LOG("getMeshVertexSize=" << p.getMeshVertexSize(s.myMesh));

  bool initialData = p.requiresInitialData();
  LOG("requiresInitialData=" << initialData);
  if (initialData) {
    p.writeData(s.myMesh, s.writeData, ids, makeValues(s, -1));
    LOG("wrote initial data");
  }

  p.initialize();
  LOG("initialize done");
  LOG("isCouplingOngoing=" << p.isCouplingOngoing());
  LOG("isTimeWindowComplete=" << p.isTimeWindowComplete());
  LOG("getMaxTimeStepSize=" << p.getMaxTimeStepSize());

  int    step   = 0;
  int    window = 0;
  int    guard  = 0;
  double subDt  = -1.0;
  while (p.isCouplingOngoing()) {
    if (++guard > 100) {
      LOG("GUARD tripped: too many steps");
      return 1;
    }
    bool wc = p.requiresWritingCheckpoint();
    LOG("requiresWritingCheckpoint=" << wc);

    double maxDt = p.getMaxTimeStepSize();
    LOG("loop maxDt=" << maxDt);
    // subcycle: take two half-window steps in the second time window
    double dt = maxDt;
    if (fixedDt > 0.0) {
      dt = std::min(fixedDt, maxDt);
    } else if (subcycle && window == 1) {
      if (subDt < 0.0)
        subDt = maxDt * 0.5; // fixed substep, set when entering the window
      dt = std::min(subDt, maxDt);
    }

    std::vector<double> rv(3 * s.readDims);
    p.readData(s.myMesh, s.readData, ids, 0.0, rv);
    LOG("V read@0=" << fmtv(rv));
    p.readData(s.myMesh, s.readData, ids, dt, rv);
    LOG("V read@dt=" << fmtv(rv));

    auto wv = makeValues(s, step);
    p.writeData(s.myMesh, s.writeData, ids, wv);
    LOG("V wrote=" << fmtv(wv));

    p.advance(dt);
    LOG("advance(" << (dt == maxDt ? "maxDt" : "half") << ") ok");

    bool rc = p.requiresReadingCheckpoint();
    LOG("requiresReadingCheckpoint=" << rc);
    bool twc = p.isTimeWindowComplete();
    LOG("isTimeWindowComplete=" << twc);
    if (twc) {
      ++window;
      subDt = -1.0;
    }
    ++step;
  }
  LOG("loop done, steps=" << step << " windows=" << window);
  LOG("final getMaxTimeStepSize=" << p.getMaxTimeStepSize());
  std::vector<double> rv(3 * s.readDims);
  p.readData(s.myMesh, s.readData, ids, 0.0, rv);
  LOG("V final read@0=" << fmtv(rv));
  p.finalize();
  LOG("finalize done");
  return 0;
}

// ---------------------------------------------------------------------------
// Error scenarios. "preinit", "ctor-badrank" and "emptymesh" don't need a partner.
static int runErrors(const std::string &name, const std::string &config, const std::string &which)
{
  Setup s = makeSetup(name);

  if (which == "ctor-badrank") {
    expectError("ctor-badrank", [&] { Participant p(name, config, 2, 1); });
    expectError("ctor-badsize", [&] { Participant p(name, config, 0, 0); });
    expectError("ctor-negrank", [&] { Participant p(name, config, -1, 2); });
    expectError("ctor-nullcomm", [&] { Participant p(name, config, 0, 1, nullptr); });
    expectError("ctor-emptyname", [&] { Participant p("", config, 0, 1); });
    expectError("ctor-unknown-participant", [&] { Participant p("NoSuchSolver", config, 0, 1); });
    return 0;
  }

  if (which == "preinit") {
    Participant p(name, config, 0, 1);
    expectError("getMeshDimensions-unknown", [&] { p.getMeshDimensions("NopeMesh"); });
    expectError("getDataDimensions-unknown", [&] { p.getDataDimensions(s.myMesh, "NopeData"); });
    expectError("getDataDimensions-unknownmesh", [&] { p.getDataDimensions("NopeMesh", s.writeData); });
    LOG("getMeshVertexSize-unknown-result: " << [&]() -> std::string {
      try { return std::to_string(p.getMeshVertexSize("NopeMesh")); }
      catch (const std::exception &e) { return std::string("EXC ") + e.what(); }
    }());
    LOG("requiresMeshConnectivityFor-unknown-result: " << [&]() -> std::string {
      try { return std::to_string(p.requiresMeshConnectivityFor("NopeMesh")); }
      catch (const std::exception &e) { return std::string("EXC ") + e.what(); }
    }());
    LOG("getMeshVertexSize-received-preinit: " << [&]() -> std::string {
      try { return std::to_string(p.getMeshVertexSize(s.otherMesh)); }
      catch (const std::exception &e) { return std::string("EXC ") + e.what(); }
    }());
    expectError("setMeshVertex-wrongdims", [&] {
      std::vector<double> pos = {1.0, 2.0, 3.0}; // mesh is 2D
      p.setMeshVertex(s.myMesh, pos);
    });
    expectError("setMeshVertices-mismatch", [&] {
      std::vector<double>   c = {0.0, 0.0, 1.0}; // 3 coords for 2 ids in 2D
      std::vector<VertexID> ids(2);
      p.setMeshVertices(s.myMesh, c, ids);
    });
    expectError("readData-beforeinit", [&] {
      std::vector<VertexID> ids = {0};
      std::vector<double>   v(1 * s.readDims);
      p.readData(s.myMesh, s.readData, ids, 0.0, v);
    });
    expectError("isCouplingOngoing-beforeinit", [&] { p.isCouplingOngoing(); });
    expectError("getMaxTimeStepSize-beforeinit", [&] { p.getMaxTimeStepSize(); });
    expectError("isTimeWindowComplete-beforeinit", [&] { p.isTimeWindowComplete(); });
    expectError("requiresWritingCheckpoint-beforeinit", [&] { p.requiresWritingCheckpoint(); });
    expectError("requiresReadingCheckpoint-beforeinit", [&] { p.requiresReadingCheckpoint(); });
    expectError("stopProfiling-unbalanced", [&] { p.stopLastProfilingSection(); });
    expectError("startProfiling-slash", [&] { p.startProfilingSection("a/b"); });

    // now set up mesh, then write-related errors (writeData is legal pre-init)
    auto ids = setupMesh(p, s);
    expectError("writeData-notwritable", [&] {
      std::vector<double> v(3 * s.readDims, 1.0);
      p.writeData(s.myMesh, s.readData, ids, v);
    });
    expectError("writeData-unknowndata", [&] {
      std::vector<double> v(3, 1.0);
      p.writeData(s.myMesh, "NopeData", ids, v);
    });
    expectError("writeData-wrongsize", [&] {
      std::vector<double> v(3 * s.writeDims + 1, 1.0);
      p.writeData(s.myMesh, s.writeData, ids, v);
    });
    expectError("writeData-invalidid", [&] {
      std::vector<VertexID> bad = {0, 1, 99};
      std::vector<double>   v(3 * s.writeDims, 1.0);
      p.writeData(s.myMesh, s.writeData, bad, v);
    });
    expectError("setMeshEdge-invalidid", [&] { p.setMeshEdge(s.myMesh, 0, 99); });
    expectError("setMeshTriangle-repeated", [&] { p.setMeshTriangle(s.myMesh, 0, 1, 1); });
    expectError("setMeshEdges-odd", [&] {
      std::vector<VertexID> e = {0, 1, 2};
      p.setMeshEdges(s.myMesh, e);
    });
    expectError("resetMesh-not-enabled", [&] { p.resetMesh(s.myMesh); });
    expectError("setMeshAccessRegion-not-received", [&] {
      std::vector<double> bb = {0.0, 1.0, 0.0, 1.0};
      p.setMeshAccessRegion(s.myMesh, bb);
    });
    expectError("writeGradientData-nogradient", [&] {
      std::vector<double> g(3 * 2 * s.writeDims, 0.0);
      p.writeGradientData(s.myMesh, s.writeData, ids, g);
    });
    expectError("setMeshVertex-unknownmesh", [&] {
      std::vector<double> pos = {0.5, 0.5};
      p.setMeshVertex("NopeMesh", pos);
    });
    expectError("setMeshVertices-unknownmesh", [&] {
      std::vector<double>   c = {0.0, 0.0, 1.0, 1.0};
      std::vector<VertexID> vids(2);
      p.setMeshVertices("NopeMesh", c, vids);
    });
    expectError("setMeshEdge-unknownmesh", [&] { p.setMeshEdge("NopeMesh", 0, 1); });
    // last: this aborts (PRECICE_ASSERT) in assertion-enabled builds of the real
    // library, so nothing may follow
    expectError("advance-beforeinit", [&] { p.advance(0.1); });
    return 0;
  }

  if (which == "emptymesh") {
    Participant p(name, config, 0, 1);
    expectError("initialize-emptymesh", [&] { p.initialize(); });
    return 0;
  }

  if (which == "postinit") {
    Participant p(name, config, 0, 1);
    auto        ids = setupMesh(p, s);
    p.initialize();
    LOG("initialize done");
    expectError("initialize-twice", [&] { p.initialize(); });
    expectError("requiresInitialData-afterinit", [&] { p.requiresInitialData(); });
    expectError("setMeshVertex-afterinit", [&] {
      std::vector<double> pos = {0.5, 0.5};
      p.setMeshVertex(s.myMesh, pos);
    });
    double maxDt = p.getMaxTimeStepSize();
    expectError("readData-negative-time", [&] {
      std::vector<double> v(3 * s.readDims);
      p.readData(s.myMesh, s.readData, ids, -0.1, v);
    });
    expectError("readData-beyond-window", [&] {
      std::vector<double> v(3 * s.readDims);
      p.readData(s.myMesh, s.readData, ids, maxDt * 1.5, v);
    });
    expectError("readData-notreadable", [&] {
      std::vector<double> v(3 * s.writeDims);
      p.readData(s.myMesh, s.writeData, ids, 0.0, v);
    });
    expectError("readData-wrongsize", [&] {
      std::vector<double> v(3 * s.readDims + 1);
      p.readData(s.myMesh, s.readData, ids, 0.0, v);
    });
    expectError("readData-invalidid", [&] {
      std::vector<VertexID> bad = {0, 1, 99};
      std::vector<double>   v(3 * s.readDims);
      p.readData(s.myMesh, s.readData, bad, 0.0, v);
    });
    // valid step so the partner (running the same scenario) can proceed
    std::vector<double> wv = makeValues(s, 0);
    p.writeData(s.myMesh, s.writeData, ids, wv);
    p.advance(maxDt);
    LOG("advance ok after error probes");
    p.finalize();
    LOG("finalize done");
    expectError("finalize-twice", [&] { p.finalize(); });
    expectError("isCouplingOngoing-afterfinalize", [&] { p.isCouplingOngoing(); });
    expectError("readData-afterfinalize", [&] {
      std::vector<double> v(3 * s.readDims);
      p.readData(s.myMesh, s.readData, ids, 0.0, v);
    });
    // last: may abort in assertion-enabled builds of the real library
    expectError("advance-afterfinalize", [&] { p.advance(0.1); });
    return 0;
  }

  if (which == "afterloop") {
    // run coupling to completion, then check post-coupling behavior
    Participant p(name, config, 0, 1);
    auto        ids = setupMesh(p, s);
    p.initialize();
    int step = 0, guard = 0;
    while (p.isCouplingOngoing()) {
      if (++guard > 100)
        return 1;
      double              maxDt = p.getMaxTimeStepSize();
      std::vector<double> wv    = makeValues(s, step);
      p.writeData(s.myMesh, s.writeData, ids, wv);
      p.advance(maxDt);
      ++step;
    }
    LOG("coupling done after " << step << " steps");
    expectError("advance-after-coupling-done", [&] { p.advance(0.1); });
    expectError("writeData-after-coupling-done", [&] {
      std::vector<double> wv = makeValues(s, step);
      p.writeData(s.myMesh, s.writeData, ids, wv);
    });
    expectError("readData-nonzero-after-done", [&] {
      std::vector<double> v(3 * s.readDims);
      p.readData(s.myMesh, s.readData, ids, 0.5, v);
    });
    std::vector<double> v(3 * s.readDims);
    p.readData(s.myMesh, s.readData, ids, 0.0, v);
    LOG("read@0 after done ok");
    LOG("getMaxTimeStepSize after done=" << p.getMaxTimeStepSize());
    p.finalize();
    LOG("finalize done");
    return 0;
  }

  // adapter "forgets" the checkpoint API in implicit coupling: real preCICE
  // must reject the first advance(); run with the implicit config
  if (which == "nocheckpoint") {
    Participant p(name, config, 0, 1);
    auto        ids = setupMesh(p, s);
    if (p.requiresInitialData()) {
      p.writeData(s.myMesh, s.writeData, ids, makeValues(s, -1));
    }
    p.initialize();
    LOG("initialize done");
    double              maxDt = p.getMaxTimeStepSize();
    std::vector<double> wv    = makeValues(s, 0);
    p.writeData(s.myMesh, s.writeData, ids, wv);
    expectError("advance-without-checkpoint", [&] { p.advance(maxDt); });
    std::cout.flush();
    std::_Exit(0);
  }

  // single advance probes: real preCICE cannot recover from a failed advance,
  // so each of these runs in a fresh process (partner runs lifecycle and dies with us)
  if (which == "advzero" || which == "advneg" || which == "advinf" || which == "advbig") {
    Participant p(name, config, 0, 1);
    auto        ids = setupMesh(p, s);
    (void) ids;
    p.initialize();
    LOG("initialize done");
    double maxDt = p.getMaxTimeStepSize();
    if (which == "advzero")
      expectError("advance-zero", [&] { p.advance(0.0); });
    else if (which == "advneg")
      expectError("advance-negative", [&] { p.advance(-0.5); });
    else if (which == "advinf")
      expectError("advance-inf", [&] { p.advance(INFINITY); });
    else
      expectError("advance-toolarge", [&] { p.advance(maxDt * 1.5); });
    std::cout.flush();
    std::_Exit(0); // skip destructor: state is unrecoverable in the real library
  }

  std::cerr << "unknown error scenario " << which << std::endl;
  return 2;
}

int main(int argc, char **argv)
{
  if (argc < 4) {
    std::cerr << "usage: driver <participant> <config> <scenario> [sub-scenario]" << std::endl;
    return 2;
  }
  std::string name     = argv[1];
  std::string config   = argv[2];
  std::string scenario = argv[3];

  try {
    if (scenario == "lifecycle")
      return runLifecycle(name, config, /*subcycle=*/false);
    if (scenario == "lifecycle-sub")
      return runLifecycle(name, config, /*subcycle=*/true);
    if (scenario == "lifecycle-fixed")
      return runLifecycle(name, config, /*subcycle=*/false, /*fixedDt=*/0.4);
    if (scenario == "errors")
      return runErrors(name, config, argc > 4 ? argv[4] : "preinit");
  } catch (const std::exception &e) {
    LOG("FATAL " << e.what());
    return 1;
  }
  std::cerr << "unknown scenario " << scenario << std::endl;
  return 2;
}
