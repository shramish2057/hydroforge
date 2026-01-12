# Input Data

This page describes the input data formats accepted by HydroForge.

## Scenario Configuration (TOML)

The scenario file is a TOML configuration that defines all simulation inputs.

### Complete Example

```toml
[scenario]
name = "Urban Flood Simulation"
description = "100-year storm event simulation"

[grid]
nx = 256                # Cells in x-direction
ny = 256                # Cells in y-direction
dx = 5.0                # Cell width (m)
dy = 5.0                # Cell height (m), defaults to dx
x_origin = 0.0          # X coordinate of lower-left corner
y_origin = 0.0          # Y coordinate of lower-left corner
crs = "EPSG:32632"      # Coordinate reference system

[input]
dem = "elevation.txt"   # DEM file path (relative to scenario.toml)
roughness = "manning.txt"  # Roughness file OR scalar value
rainfall = "storm.csv"  # Rainfall time series

[parameters]
t_end = 7200.0          # Simulation duration (s)
dt_max = 2.0            # Maximum timestep (s)
cfl = 0.7               # CFL number (0.5-0.9)
h_min = 0.001           # Minimum wet depth (m)
output_interval = 300.0 # Output interval (s)

[output]
directory = "results"   # Output directory
points = [[128, 128], [64, 64], [192, 192]]  # Monitor points [i, j]
```

## Digital Elevation Model (DEM)

### Text Format (.txt)

Space-separated matrix of elevation values:

```
12.5 12.4 12.3 12.2 12.1
12.4 12.3 12.2 12.1 12.0
12.3 12.2 12.1 12.0 11.9
```

- First row corresponds to j=1 (bottom of domain)
- First column corresponds to i=1 (left of domain)
- Values in meters above datum

### GeoTIFF Format (.tif)

Requires ArchGDAL.jl. Standard single-band GeoTIFF with:

- Float32 or Float64 values
- Proper georeferencing
- NoData handling

### Requirements

- Dimensions must match grid (nx × ny)
- No NaN or Inf values (except NoData cells)
- Units: meters above datum

## Roughness Data

### Scalar Value

For uniform roughness:
```toml
roughness = 0.03  # Manning's n everywhere
```

### Spatially Variable

Text file with same dimensions as DEM:
```toml
roughness = "manning.txt"
```

### Value Ranges

- Minimum: 0.01 (very smooth)
- Maximum: 0.50 (heavily vegetated)
- Typical urban: 0.02-0.05

## Rainfall Data

### CSV Format

```csv
time,intensity
0,0
300,5.0
600,15.0
900,30.0
1200,50.0
1500,30.0
1800,10.0
2100,0
```

- **time**: Seconds from simulation start
- **intensity**: mm/hr (millimeters per hour)
- Linear interpolation between points
- Must be monotonically increasing in time

### Design Storm Examples

**1-hour intense storm (50 mm total):**
```csv
time,intensity
0,0
600,20
1200,50
1800,50
2400,30
3000,10
3600,0
```

**6-hour moderate storm:**
```csv
time,intensity
0,0
3600,15
7200,25
10800,20
14400,15
18000,5
21600,0
```

## Validation

HydroForge validates inputs automatically:

```julia
using HydroForge

# Validate a scenario
scenario = load_scenario_from_toml("scenario.toml")
validate_scenario(scenario)  # Throws on invalid data, warns on unusual values
```

### Common Validation Errors

| Error | Cause | Solution |
|-------|-------|----------|
| DimensionMismatch | DEM size ≠ grid size | Check nx, ny match DEM |
| ArgumentError: negative Manning | n ≤ 0 | Use positive values |
| NaN in DEM | Missing data | Replace with valid elevations |
