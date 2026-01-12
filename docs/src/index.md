# HydroForge.jl

*Real-Time Urban Flood Risk Simulator*

HydroForge is a high-performance, open-source Julia package for simulating urban stormwater runoff and surface inundation using the shallow water equations.

## Features

- **Real-time simulation**: Efficient local inertial solver for interactive analysis
- **Urban focus**: Designed for pluvial flood prediction in cities
- **Easy to use**: Simple API for running simulations
- **Extensible**: Modular architecture for customization
- **Well-tested**: Comprehensive test suite with analytical benchmarks

## Quick Example

```julia
using HydroForge

# Run the bundled demo scenario
HydroForge.run_demo()

# Or run your own scenario
HydroForge.run("path/to/scenario.toml")
```

## Installation

```julia
using Pkg
Pkg.add("HydroForge")
```

See [Installation](@ref) for detailed instructions.

## Documentation

```@contents
Pages = [
    "installation.md",
    "quickstart.md",
    "concepts.md",
]
Depth = 1
```

## Mathematical Foundation

HydroForge uses the **local inertial approximation** of the shallow water equations:

**Continuity:**
```math
\frac{\partial h}{\partial t} + \frac{\partial q_x}{\partial x} + \frac{\partial q_y}{\partial y} = R
```

**Momentum (x-direction):**
```math
\frac{\partial q_x}{\partial t} + gh\frac{\partial \eta}{\partial x} + \frac{gn^2|q_x|q_x}{h^{10/3}} = 0
```

Where:
- ``h`` = water depth (m)
- ``q`` = unit discharge (m²/s)
- ``\eta`` = water surface elevation (m)
- ``R`` = rainfall rate (m/s)
- ``n`` = Manning's roughness coefficient
- ``g`` = gravitational acceleration (m/s²)

## License

HydroForge is released under the [MIT License](https://opensource.org/licenses/MIT).

## Acknowledgments

HydroForge builds upon ideas from:
- [ShallowWaters.jl](https://github.com/milankl/ShallowWaters.jl)
- [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)
- The broader Julia scientific computing community
