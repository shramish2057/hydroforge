# Quick Start

This guide will get you running your first flood simulation in 5 minutes.

## Running the Demo

The simplest way to start is with the bundled demo:

```julia
using HydroForge

# Run the demo scenario
HydroForge.run_demo()
```

This simulates a 1-hour rainfall event on a 64×64 synthetic urban catchment.

### Demo Output

The demo creates results in `assets/demo_data/runs/<run_id>/`:

| File | Description |
|------|-------------|
| `max_depth.txt` | Maximum water depth at each cell (m) |
| `arrival_time.txt` | Time when water first arrived (s) |
| `max_velocity.txt` | Maximum flow velocity (m/s) |
| `hydrograph_*.csv` | Time series at output points |
| `metadata.json` | Run metadata and statistics |

## Creating Your Own Scenario

### Step 1: Prepare Input Data

You need three inputs:

1. **DEM** (Digital Elevation Model): Grid of ground elevations
2. **Roughness**: Manning's n values for each cell
3. **Rainfall**: Time series of rainfall intensity

### Step 2: Create Scenario File

Create a `scenario.toml` file:

```toml
[scenario]
name = "My First Scenario"
description = "Custom urban flood simulation"

[grid]
nx = 100          # Number of cells in x
ny = 100          # Number of cells in y
dx = 10.0         # Cell size in meters

[input]
dem = "dem.txt"           # Path to elevation data
roughness = 0.03          # Uniform roughness (or path to file)
rainfall = "rainfall.csv" # Path to rainfall data

[parameters]
t_end = 3600.0            # Simulation duration (seconds)
dt_max = 1.0              # Maximum timestep
cfl = 0.7                 # CFL number (stability)
output_interval = 60.0    # Output every 60 seconds

[output]
directory = "results"
points = [[50, 50], [25, 75]]  # Monitor points (i, j)
```

### Step 3: Prepare DEM

Create `dem.txt` with space-separated elevation values:

```
10.5 10.4 10.3 10.2 ...
10.4 10.3 10.2 10.1 ...
...
```

### Step 4: Prepare Rainfall

Create `rainfall.csv`:

```csv
time,intensity
0,0
600,10
1200,30
1800,50
2400,30
3000,10
3600,0
```

Time is in seconds, intensity in mm/hr.

### Step 5: Run Simulation

```julia
using HydroForge

results = HydroForge.run("path/to/scenario.toml")

println("Simulation completed!")
println("Maximum depth: ", results["max_depth"], " m")
println("Results saved to: ", results["output_dir"])
```

## Visualizing Results

```julia
using DelimitedFiles
using Plots

# Load max depth
max_depth = readdlm("results/max_depth.txt")

# Create heatmap
heatmap(max_depth',
    title="Maximum Flood Depth",
    xlabel="X (cells)",
    ylabel="Y (cells)",
    colorbar_title="Depth (m)",
    c=:blues)

savefig("flood_map.png")
```

## Next Steps

- [Concepts](@ref): Understand the physics and numerics
- [Input Data](@ref): Detailed input data specifications
- [API Reference](@ref): Full API documentation
