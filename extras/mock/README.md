# preCICE Mock Participant

A lightweight mock implementation of the preCICE Participant API for testing and development.

## Features

- **Full API Coverage**: Implements the complete preCICE Participant API
- **Configuration Validation**: Parses preCICE XML configs and validates API calls
- **Implicit Coupling Support**: Handles checkpoint requirements for implicit coupling schemes
- **Flexible Data Exchange**: Three modes for handling read/write data operations
- **Error Handling**: Throws proper preCICE exceptions with helpful error messages

## Data Exchange Modes

The mock participant supports three modes for `readData()` operations, configured via an optional `mock-config.xml` file:

### 1. Random Mode
Returns random seeded data (useful for testing error handling).

```xml
<data-mapping mesh="MeshName" data="DataName" mode="random" />
```

### 2. Buffer Mode (Default)
Returns the data previously written via `writeData()`. This mimics the behavior of the Python adapter mock and is the default if no mock-config is provided.

```xml
<data-mapping mesh="MeshName" data="DataName" mode="buffer" />
```

### 3. Scaled Buffer Mode
Returns the buffered write data multiplied by a scalar or element-wise by a vector.

**Scalar multiplication:**
```xml
<data-mapping mesh="MeshName" data="DataName" mode="scaled">
  <scalar-multiplier value="2.0" />
</data-mapping>
```

**Element-wise vector multiplication:**
```xml
<data-mapping mesh="MeshName" data="DataName" mode="scaled">
  <vector-multiplier values="1.0;2.0;3.0" />
</data-mapping>
```

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
<mock-config>
  <!-- Random data for temperature -->
  <data-mapping mesh="FluidMesh" data="Temperature" mode="random" />

  <!-- Return written data as-is for pressure -->
  <data-mapping mesh="FluidMesh" data="Pressure" mode="buffer" />

  <!-- Scale displacement by factor of 2 -->
  <data-mapping mesh="StructureMesh" data="Displacement" mode="scaled">
    <scalar-multiplier value="2.0" />
  </data-mapping>

  <!-- Scale force vector components differently -->
  <data-mapping mesh="StructureMesh" data="Force" mode="scaled">
    <vector-multiplier values="1.0;1.5;2.0" />
  </data-mapping>
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

The mock is built as part of the preCICE library. Link against `libpreciceMocked.so` instead of `libprecice.so`:

```cmake
target_link_libraries(your_solver preciceMocked)
```

Alternatively build normally and use LD_PRELOAD to load the mock at runtime

```
LD_PRELOAD=/path/to/precicedir/precice/build/libpreciceMocked.so.3.3.0
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
2. Optionally create a mock-config.xml for specific data behavior
3. Run your adapter against the mock
4. Verify correct API usage patterns
5. Test error handling by introducing configuration mistakes
6. Validate checkpoint handling for implicit coupling

The mock is ideal for:
- Unit testing adapters without running full simulations
- Testing error handling code paths
- Developing adapters without a coupling partner
- CI/CD integration testing
