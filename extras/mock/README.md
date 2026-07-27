# preCICE Mock Participant

A lightweight mock implementation of the preCICE Participant API for testing and developing adapters: run your solver against the mock — without a coupling partner and without the full preCICE stack — to validate API usage, exercise error paths, and run adapter tests in CI.

## Building

The mock is built as part of the normal preCICE build and produces the shared library `libpreciceMocked.so`.

To build only the mock, but not `libprecice`, `testprecice`, etc.:

```bash
cmake --build . --target preciceMocked   # or: make preciceMocked
```

### Dependencies

The mock needs a proper subset of the [dependencies of preCICE](https://precice.org/installation-source-dependencies.html):

- CMake and a C++17 compiler
- Eigen
- libxml2
- MPI (optional, following the `PRECICE_FEATURE_MPI_COMMUNICATION` option)

Boost, PETSc, and network communication libraries are not needed.

The mock does not link against `libprecice`; it re-implements the `Participant` API in [preciceMocked.cpp](preciceMocked.cpp). From preCICE itself it only uses a few headers (exception types, type aliases, the `span` shim, version constants) and the bundled fmt library.

## How to use

The mock reads your standard preCICE configuration file. The configuration must be valid — the mock assumes a configuration that passes `precice-tools check` and reports undefined behavior on invalid ones.

There are two ways to run a solver against the mock:

### Option 1: LD_PRELOAD (no rebuild needed)

Build your solver against the real preCICE as usual, then preload the mock at runtime. This also works for solvers using the language bindings, since these load `libprecice` underneath.

Running a Python solver (e.g. a FEniCS case using the FEniCS adapter):

```bash
LD_PRELOAD=/path/to/precice/build/libpreciceMocked.so python3 solver.py
```

Running an OpenFOAM case via its run script:

```bash
LD_PRELOAD=/path/to/precice/build/libpreciceMocked.so ./run.sh
```

### Option 2: Link against the mock

Link `libpreciceMocked.so` instead of `libprecice.so`:

```cmake
target_link_libraries(your_solver preciceMocked)
```

### Testing workflow

1. Run your adapter with your preCICE configuration against the mock
2. Verify correct API usage patterns
3. Test error handling by introducing configuration mistakes
4. Validate checkpoint handling for implicit coupling
5. Optionally, add a `precice-mock-config.xml` to control the returned data (see [Configuration](#configuration))

### Usage example

```cpp
#include <precice/precice.hpp>

int main() {
  precice::Participant participant("SolverOne", "precice-config.xml", 0, 1);

  // Initialize reads precice-config.xml and precice-mock-config.xml (if present)
  participant.initialize();

  // Set up mesh
  std::vector<double> coords = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
  std::vector<VertexID> ids(2);
  participant.setMeshVertices("MyMesh", coords, ids);

  // Write some data
  std::vector<double> writeValues = {1.0, 2.0};
  participant.writeData("MyMesh", "Temperature", ids, writeValues);

  participant.advance(0.1);

  // Read data - behavior depends on precice-mock-config.xml
  std::vector<double> readValues(2);
  participant.readData("MyMesh", "Temperature", ids, 0.0, readValues);
  // readValues will contain:
  // - Random data if mode="random"
  // - {1.0, 2.0} if mode="buffer" (default)
  // - {2.0, 4.0} if mode="scaled" with scalar=2.0

  participant.finalize();
  return 0;
}
```

## Behavior

- **Full API coverage**: Implements the complete preCICE Participant API
- **Caller validation**: Validates every API call against the preCICE configuration — participant, mesh, and data names must exist, read/write directions must match, and steering calls must arrive in a legal order
- **Data exchange**: `readData()` returns buffered, scaled, or random data depending on the mock configuration (see [Configuration](#configuration)); `relativeReadTime` is ignored
- **Implicit coupling**: Checkpoint requirements are enforced; convergence measures are not evaluated — the scheme converges after a fixed number of iterations (default 2, see [Iteration override](#iteration-override))
- **Termination**: The mock terminates based on `max-time` / `max-time-windows` from the configuration. If the configuration defines neither (which is valid — real preCICE would run until stopped by the coupling partner), the mock has no partner to stop it, so it adds `max-time-windows="100"` and prints a warning. With `<time-window-size method="first-participant" />`, every participant prescribes the window size, since the mock has no partner that could prescribe it

### Error handling

The mock validates all API calls against the configuration and throws `precice::Error` with descriptive messages:

```text
Error: Data 'Temperature' on mesh 'FluidMesh' is not configured for writing by participant 'SolverOne'.
Please add <write-data name="Temperature" mesh="FluidMesh" /> to the configuration.
```

This helps catch configuration mistakes early in development.

## Configuration

### preCICE configuration

Use your standard, valid preCICE configuration file — the mock needs no dedicated configuration.

### Mock configuration (optional)

The mock behavior can be adjusted with an optional file `precice-mock-config.xml`, placed in the same directory as the preCICE configuration file. An example can be found in this folder: [precice-mock-config.xml](precice-mock-config.xml).

#### Logging modes

The mock supports two logging modes:

- **Mock mode (default)**: Minimal output showing only essential mock behavior with simple messages prefixed by `[precice-mock]`.

  ```xml
  <logging-mode mode="mock" />
  ```

- **PrecICE mode**: Verbose output that mimics the real preCICE library, showing detailed initialization, mesh setup, and shutdown messages with the `---[precice]` prefix.

  ```xml
  <logging-mode mode="precice" />
  ```

#### Data exchange modes

The mock participant supports three modes for `readData()` operations. The just-in-time mapping API participates in the same machinery: `mapAndReadData()` follows the same modes as `readData()`, and `writeAndMapData()` fills the same per mesh/data buffers as `writeData()` (but is exempt from the read-size check, since just-in-time calls may use a different vertex count on every call).

**1. Buffer mode (default)**: Returns the data previously written via `writeData()`. This is the default if no mock-config is provided.

Buffers are kept per mesh/data pair. A read looks up the buffer for the same mesh and data first, then falls back to the same data name written on another mesh, then to the last written data of any kind (values repeat cyclically if the buffer is shorter than the read).

```xml
<mocked-data mesh="MeshName" data="DataName" mode="buffer" />
```

**2. Random mode**: Returns random seeded data (useful for testing error handling). Optional bounds and seed can be specified with nested elements (bounds defaults: 0.0 to 1.0, seed defaults to rank-based).

```xml
<mocked-data mesh="MeshName" data="DataName" mode="random">
  <bounds lower="-1.0" upper="1.0" />
  <seed value="42" />
</mocked-data>
```

**3. Scaled buffer mode**: Returns the buffered write data multiplied by a scalar or element-wise by a vector.

Scalar multiplication:

```xml
<mocked-data mesh="MeshName" data="DataName" mode="scaled">
  <scalar-multiplier value="2.0" />
</mocked-data>
```

Element-wise vector multiplication:

```xml
<mocked-data mesh="MeshName" data="DataName" mode="scaled">
  <vector-multiplier values="1.0;2.0;3.0" />
</mocked-data>
```

The vector multiplier is applied cyclically when it is shorter than the data arity, so you can provide one value per component or a shorter repeating pattern.

A global default mode and multipliers can be set for all data items not explicitly configured:

```xml
<mock-config>
  <!-- Set default mode for all unconfigured data -->
  <mocked-data-default mode="scaled">
    <scalar-multiplier value="1.5" />
  </mocked-data-default>

  <!-- Specific configs override the global default -->
  <mocked-data mesh="MeshOne" data="SpecialData" mode="random">
    <bounds lower="0.0" upper="5.0" />
    <seed value="123" />
  </mocked-data>
</mock-config>
```

In this example, `SpecialData` uses random mode with custom bounds and seed, while all other data items use scaled mode with multiplier 1.5. If no `<mocked-data-default>` is specified, the global default is buffer mode.

#### Iteration override

For implicit coupling, you can override the `max-iterations` value from the preCICE configuration to converge faster during testing, without modifying the preCICE configuration file:

```xml
<mock-config>
  <!-- Value: positive integer for custom iterations, -1 to use config value -->
  <max-iterations-override value="2" />
</mock-config>
```

If you omit the tag entirely, the mock defaults to `2` iterations. If you set the value to `-1`, the mock uses the `max-iterations` value from the preCICE configuration. If the override value is higher than the configured value, a warning is printed but the override is still respected.

## Limitations

- **No real communication**: Each rank mocks the coupling independently — data read back comes from the rank-local buffers. An MPI communicator passed to the constructor is accepted (and validated), but not used for communication. Consequently, `getMeshVertexSize()` on a received mesh without API access returns 0.
- **Gradient data is discarded**: `writeGradientData()` validates its arguments, but the gradients are not stored — there is no read counterpart in the preCICE API (yet) that could return them.
- **Connectivity is not stored**: `setMeshEdges()`, `setMeshTriangles()`, and related calls validate their arguments but are otherwise no-ops.
- **Multiple coupling schemes are merged**: With several coupling schemes in one configuration, the schemes this participant belongs to are merged into one (largest termination bounds win).
