# Installation

## Requirements

- **Julia**: Version 1.10 or higher
- **Operating System**: Linux, macOS, or Windows
- **Memory**: Minimum 4GB RAM (8GB+ recommended for large simulations)

## Installing Julia

If you don't have Julia installed:

1. Download from [julialang.org](https://julialang.org/downloads/)
2. Follow the installation instructions for your platform
3. Verify installation: `julia --version`

## Installing HydroForge

### From Julia Package Registry (Recommended)

```julia
using Pkg
Pkg.add("HydroForge")
```

### From Source (Development)

```julia
using Pkg
Pkg.develop(url="https://github.com/hydroforge/HydroForge.git")
```

### Verifying Installation

```julia
using HydroForge

# Check version
HydroForge.version()

# Display package info
HydroForge.info()

# Run demo to verify everything works
HydroForge.run_demo()
```

## Optional Dependencies

### For GeoTIFF Support

To read/write GeoTIFF files (recommended for real-world applications):

```julia
using Pkg
Pkg.add("ArchGDAL")
```

Note: ArchGDAL requires GDAL libraries. See [ArchGDAL documentation](https://yeesian.com/ArchGDAL.jl/stable/) for platform-specific instructions.

### For Visualization

```julia
using Pkg
Pkg.add("Plots")
Pkg.add("GeoMakie")  # For geospatial visualization
```

## Troubleshooting

### Common Issues

**Package not found:**
```julia
# Update registry
Pkg.Registry.update()
Pkg.add("HydroForge")
```

**Dependency conflicts:**
```julia
# Create a fresh environment
Pkg.activate("my_project")
Pkg.add("HydroForge")
```

**Demo fails to run:**
- Ensure you have write permissions in the current directory
- Check that all dependencies are installed: `Pkg.instantiate()`

### Getting Help

- Check [GitHub Issues](https://github.com/hydroforge/HydroForge/issues)
- Open a new issue with your error message and environment details
