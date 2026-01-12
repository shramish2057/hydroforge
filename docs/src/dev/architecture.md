# Architecture

## Overview

HydroForge follows a modular architecture with clear separation of concerns:

```
HydroForge/
├── src/
│   ├── HydroForge.jl      # Main module
│   ├── config/            # Configuration & defaults
│   ├── types/             # Core data types
│   ├── numerics/          # Numerical methods
│   ├── physics/           # Physical processes
│   ├── io/                # Input/output
│   ├── models/            # Simulation models
│   ├── cli/               # Command-line interface
│   └── api/               # HTTP API server
├── test/                  # Test suite
├── docs/                  # Documentation
├── benchmark/             # Performance benchmarks
└── assets/                # Demo data
```

## Module Structure

### Types (`src/types/`)

Core data structures:

- **`grid.jl`**: Computational grid definition
- **`topography.jl`**: Terrain and roughness data
- **`state.jl`**: Simulation state (h, qx, qy)
- **`parameters.jl`**: Runtime parameters
- **`scenario.jl`**: Complete scenario bundle

### Numerics (`src/numerics/`)

Numerical algorithms:

- **`timestep.jl`**: CFL-based timestep control
- **`flux.jl`**: Local inertial flux computation
- **`boundaries.jl`**: Boundary condition handling

### Physics (`src/physics/`)

Physical process models:

- **`friction.jl`**: Manning friction
- **`rainfall.jl`**: Rainfall source term
- **`infiltration.jl`**: Infiltration (planned)

### IO (`src/io/`)

Data handling:

- **`readers.jl`**: Input file parsing
- **`writers.jl`**: Output file generation
- **`validation.jl`**: Input validation

### Models (`src/models/`)

High-level simulation control:

- **`surface2d.jl`**: 2D surface flow solver
- **`simulation.jl`**: Simulation runner

## Data Flow

```
Input Files → Readers → Scenario → Solver → Results → Writers → Output Files
                           ↓
                    SimulationState
                           ↓
                    Timestep Loop:
                    1. Compute dt (CFL)
                    2. Compute fluxes
                    3. Update depths
                    4. Apply sources
                    5. Apply boundaries
                    6. Record outputs
```

## Extension Points

### Custom Boundary Conditions

```julia
# Define new boundary type
struct MyBoundary <: BoundaryType end

# Implement apply_boundaries! for your type
function apply_boundaries!(state, ::Type{MyBoundary})
    # Custom logic
end
```

### Custom Source Terms

```julia
# Add source in simulation loop
function my_source!(h, grid, t, dt)
    # Modify h in-place
end
```

## Performance Considerations

### Memory Layout

- Arrays are column-major (Julia default)
- Access pattern: `for j in 1:ny, for i in 1:nx`
- State arrays are mutable; types are immutable

### Type Stability

- All core functions are type-stable
- Parametric types: `Grid{T}`, `SimulationState{T}`
- Type inference enables optimization

### Allocation

- Preallocated workspace: `SimulationWorkspace`
- In-place operations: `compute_flux_x!`, `update_depth!`
- Minimal allocation in hot loops
