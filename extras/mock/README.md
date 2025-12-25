# preCICE Mock Participant

A lightweight mock implementation of the preCICE Participant API for testing and development.

## Features

- **Full API Coverage**: Implements the complete preCICE Participant API
- **Caller Validation**: Parses preCICE XML configs and validates that an adapter sends API calls matching the given config
- **Implicit Coupling Support**: Handles checkpoint requirements for implicit coupling schemes
- **Flexible Data Exchange**: Three modes for handling read/write data operations
- **Error Handling**: Throws proper preCICE exceptions with helpful error messages

## Data Exchange Modes

The mock participant supports three modes for `readData()` operations, configured via an optional `mock-config.xml` file:

### 1. Random Mode
Returns random seeded data (useful for testing error handling). Optional bounds can be specified with a nested `bounds` element (defaults: 0.0 to 1.0).

```xml
<mocked-data mesh="MeshName" data="DataName" mode="random">
  <bounds lower="-1.0" upper="1.0" />
</mocked-data>
```

### 2. Buffer Mode (Default)
Returns the data previously written via `writeData()`. This is the default if no mock-config is provided.

```xml
<mocked-data mesh="MeshName" data="DataName" mode="buffer" />
```

### 3. Scaled Buffer Mode
Returns the buffered write data multiplied by a scalar or element-wise by a vector.

**Scalar multiplication:**
```xml
<mocked-data mesh="MeshName" data="DataName" mode="scaled">
  <scalar-multiplier value="2.0" />
</mocked-data>
```

**Element-wise vector multiplication:**
```xml
<mocked-data mesh="MeshName" data="DataName" mode="scaled">
  <vector-multiplier values="1.0;2.0;3.0" />
</mocked-data>
```

### Global Default Configuration

You can set a global default mode and multipliers that apply to all data items not explicitly configured:

```xml
<mock-config>
  <!-- Set default mode for all unmapped data -->
  <default mode="scaled">
    <scalar-multiplier value="1.5" />
  </default>

  <!-- Specific configs override the global default -->
  <mocked-data mesh="MeshOne" data="SpecialData" mode="random">
    <bounds lower="0.0" upper="5.0" />
  </mocked-data>
</mock-config>
```

In this example:
- `SpecialData` uses random mode
- Random mode honors optional `lower`/`upper` bounds (defaults to 0.0-1.0)
- All other data items use scaled mode with 1.5 multiplier
- If no `<default>` is specified, the global default is buffer mode

## Simulation Termination

### Explicit Coupling
For explicit coupling schemes, termination is controlled by the preCICE configuration file using `<max-time value="..."/>` or `<max-time-windows value="..."/>` in the `<coupling-scheme:.../>`.

**The mock requires at least one of these to be present to prevent infinite execution.**

### Implicit Coupling  
For implicit coupling schemes, you can configure termination in the mock-config.xml:

```xml
<termination max-implicit-rounds="5" />
```

**Parameters:**
- `max-implicit-rounds`: Number of completed implicit coupling iterations (sub-cycles) before termination. Default: 5

**Example:**
```xml
<mock-config>
  <!-- Simulate for 10 implicit rounds -->
  <termination max-implicit-rounds="10" />

  <default mode="buffer" />
  <!-- ... rest of config ... -->
</mock-config>
```

In implicit coupling, the mock counts one "implicit round" each time an entire sub-cycling iteration converges (i.e., when all iterations of a single time window are complete).

## Configuration Files

### preCICE Configuration
Use any standard preCICE configuration file. The mock will:
- Validate that the participant name exists
- Check mesh and data declarations
- Enforce read/write permissions
- Detect implicit vs explicit coupling
- Return correct dimensions for meshes and data

### Mock Configuration (Optional)
Place `mock-config.xml` in the same directory as your preCICE config file.

**Example mock-config.xml:**
```xml
<?xml version="1.0" encoding="UTF-8" ?>
<!--
  Mock Configuration for preCICE Mock Participant

  This file configures how the mock participant handles readData() calls.
  Place this file in the same directory as your precice-config.xml file.

  Three modes are supported:
  1. "random" - Returns random data (seeded)
  2. "buffer" - Returns data written by writeData() (default if no mock-config)
  3. "scaled" - Returns data written by writeData() multiplied by a scalar or vector
-->
<mock-config>

  <!-- GLOBAL DEFAULT (optional) -->
  <!-- If specified, this mode applies to all data items not explicitly configured.
       If omitted, defaults to buffer mode. -->
  <default mode="buffer">
    <!-- Optional: default scalar multiplier for scaled mode -->
    <scalar-multiplier value="1.0" />
  </default>

  <!-- Example 1: Random data mode with bounds -->
  <mocked-data mesh="FluidMesh" data="Temperature" mode="random">
    <bounds lower="0.0" upper="1.0" />
  </mocked-data>

  <!-- Example 2: Buffer mode (returns writeData values as-is) -->
  <mocked-data mesh="MeshTwo" data="Displacement" mode="buffer" />

  <!-- Example 3: Scaled buffer with scalar multiplier -->
  <mocked-data mesh="MeshOne" data="Pressure" mode="scaled">
    <scalar-multiplier value="2.0" />
  </mocked-data>

  <!-- Example 4: Scaled buffer with element-wise vector multiplier -->
  <mocked-data mesh="MeshTwo" data="Velocity" mode="scaled">
    <!-- For 3D vector data: multiply component-wise -->
    <vector-multiplier values="1.0;2.0;3.0" />
  </mocked-data>

  <!-- Example 5: Multiple scalar values can use semicolons too -->
  <mocked-data mesh="MeshThree" data="Heat-Flux" mode="scaled">
    <!-- Will cycle through multipliers for each value -->
    <vector-multiplier values="0.5;1.0;1.5;2.0" />
  </mocked-data>
</mock-config>
```

## Usage Example

```cpp
#include <precice/precice.hpp>

int main() {
  precice::Participant participant("SolverOne", "precice-config.xml", 0, 1);

  // Initialize reads precice-config.xml and mock-config.xml (if present)
  participant.initialize();

  // Set up mesh
  std::vector<double> coords = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0};
  std::vector<VertexID> ids(2);
  participant.setMeshVertices("MyMesh", coords, ids);

  // Write some data
  std::vector<double> writeValues = {1.0, 2.0};
  participant.writeData("MyMesh", "Temperature", ids, writeValues);

  participant.advance(0.1);

  // Read data - behavior depends on mock-config.xml
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

## Implicit Coupling

For implicit coupling schemes, the mock automatically:
- Returns `true` from `requiresWritingCheckpoint()` at the start of each iteration
- Returns `true` from `requiresReadingCheckpoint()` when iterating (not converged)
- Tracks iteration count and convergence status
- Reports time window completion correctly

## Building

The mock is built as part of the preCICE library.

Build normally and use LD_PRELOAD to load the mock at runtime

```
LD_PRELOAD=/path/to/precicedir/precice/build/libpreciceMocked.so.3.3.0
```

Alternatively link against `libpreciceMocked.so` instead of `libprecice.so`:

```cmake
target_link_libraries(your_solver preciceMocked)
```
## Error Handling

The mock validates all API calls against the configuration and throws `precice::Error` with descriptive messages:

```
Error: Data 'Temperature' on mesh 'FluidMesh' is not configured for writing by participant 'SolverOne'.
Please add <write-data name="Temperature" mesh="FluidMesh" /> to the configuration.
```

This helps catch configuration mistakes early in development.

## Limitations

- Does not perform actual mesh mapping or data interpolation
- Does not communicate between participants (single-process only)
- Simplified convergence logic for implicit coupling (converges after N iterations)
- No gradient data support (returns zeros)
- No mesh connectivity validation

## Testing Your Adapter

1. Create a simple preCICE config with your participant
2. Optionally, create a 'mock-config.xml' for specific data behavior
3. Run your adapter against the mock
4. Verify correct API usage patterns
5. Test error handling by introducing configuration mistakes
6. Validate checkpoint handling for implicit coupling

The mock is ideal for:
- Unit testing adapters without running full simulations
- Testing error handling code paths
- Developing adapters without a coupling partner
- CI/CD integration testing
