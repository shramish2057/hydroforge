# Output Data

## Output Directory Structure

After a simulation run, results are saved to:

```
<output_dir>/<run_id>/
├── max_depth.txt         # Maximum water depth (m)
├── arrival_time.txt      # First arrival time (s)
├── max_velocity.txt      # Maximum velocity (m/s)
├── hydrograph_i_j.csv    # Time series at point (i,j)
├── metadata.json         # Run metadata
└── run_metadata.json     # Detailed run information
```

## Raster Outputs

### Maximum Depth (`max_depth.txt`)

The maximum water depth reached at each cell during the simulation.

- **Format**: Space-separated text matrix
- **Units**: meters
- **Use**: Flood extent and depth maps

```julia
using DelimitedFiles
max_depth = readdlm("results/max_depth.txt")
```

### Arrival Time (`arrival_time.txt`)

Time when water first exceeded the threshold depth (default 1cm).

- **Format**: Space-separated text matrix
- **Units**: seconds from simulation start
- **Special value**: `Inf` for cells that never flooded

### Maximum Velocity (`max_velocity.txt`)

Maximum flow velocity magnitude at each cell.

- **Format**: Space-separated text matrix
- **Units**: m/s
- **Use**: Hazard assessment, erosion risk

## Time Series Outputs

### Hydrographs (`hydrograph_i_j.csv`)

Water depth time series at specified output points.

```csv
time,depth
0.0,0.0
60.0,0.012
120.0,0.045
...
```

- **time**: Seconds from start
- **depth**: Water depth in meters

### Loading Hydrographs

```julia
using DelimitedFiles

# Read hydrograph
data = readdlm("results/hydrograph_50_50.csv", ','; header=true)
times = data[1][:, 1]
depths = data[1][:, 2]

# Plot
using Plots
plot(times / 60, depths * 100,
    xlabel="Time (min)",
    ylabel="Depth (cm)",
    title="Hydrograph at (50, 50)")
```

## Metadata

### `metadata.json`

Summary statistics:

```json
{
    "run_id": "20240115_143052_abc123",
    "scenario_name": "Urban Flood Demo",
    "timestamp": "2024-01-15T14:30:52",
    "grid_nx": 100,
    "grid_ny": 100,
    "grid_dx": 10.0,
    "t_end": 3600.0,
    "max_depth_overall": 0.85,
    "total_rainfall_mm": 35.5
}
```

### `run_metadata.json`

Detailed run information:

```json
{
    "run_id": "20240115_143052_abc123",
    "start_time": "2024-01-15T14:30:52",
    "end_time": "2024-01-15T14:31:15",
    "julia_version": "1.10.0",
    "hydroforge_version": "0.1.0-dev",
    "git_commit": "a1b2c3d",
    "scenario_name": "Urban Flood Demo",
    "parameters": {
        "t_end": 3600.0,
        "dt_max": 1.0,
        "cfl": 0.7
    },
    "status": "completed",
    "error_message": null
}
```

## Post-Processing

### Flood Extent

```julia
# Calculate flooded area
flooded_cells = count(max_depth .> 0.1)  # Cells with >10cm
flooded_area = flooded_cells * grid.dx * grid.dy  # m²
println("Flooded area: ", flooded_area / 10000, " hectares")
```

### Volume Balance

```julia
# Check mass conservation
total_rainfall = rainfall_total * cell_area(grid) * grid.nx * grid.ny
final_volume = sum(max_depth) * cell_area(grid)
println("Rainfall volume: ", total_rainfall, " m³")
println("Final volume: ", final_volume, " m³")
```

### Visualization

```julia
using Plots

# Flood depth map
heatmap(max_depth',
    c=:blues,
    clims=(0, maximum(max_depth)),
    title="Maximum Flood Depth",
    colorbar_title="Depth (m)")

# Add contours
contour!(max_depth', levels=[0.1, 0.3, 0.5], color=:black)
```
