# HydroForge IO Readers
# Functions for reading input data files

using DelimitedFiles
using TOML

"""
    read_rainfall_csv(path::String)

Read rainfall time series from CSV file.

Expected format:
```
time,intensity
0,0
300,10.5
600,25.0
...
```

Where time is in seconds and intensity in mm/hr.

# Returns
- `RainfallEvent`: Parsed rainfall event
"""
function read_rainfall_csv(path::String)
    isfile(path) || throw(ArgumentError("File not found: $path"))

    # Read CSV data
    data = readdlm(path, ',', Float64; header=true, skipstart=0)
    matrix = data[1]

    times = matrix[:, 1]
    intensities = matrix[:, 2]

    RainfallEvent(times, intensities)
end

"""
    read_array(path::String, T=Float64)

Read a matrix from a simple text format (space-separated values).

# Arguments
- `path`: Path to file
- `T`: Element type (default Float64)

# Returns
- `Matrix{T}`: The loaded matrix
"""
function read_array(path::String, ::Type{T}=Float64) where T<:AbstractFloat
    isfile(path) || throw(ArgumentError("File not found: $path"))
    Matrix{T}(readdlm(path, T))
end

"""
    read_geotiff(path::String)

Read a GeoTIFF file and return (data, metadata).

Note: For full GeoTIFF support, install ArchGDAL.jl.
This function provides a fallback for simple array data.

# Returns
- `data::Matrix{Float64}`: Raster data
- `metadata::Dict`: Contains crs, dx, dy, x_origin, y_origin
"""
function read_geotiff(path::String)
    if endswith(path, ".tif") || endswith(path, ".tiff")
        # Check if ArchGDAL is available
        if isdefined(Main, :ArchGDAL)
            error("ArchGDAL detected but not integrated. Use ArchGDAL directly.")
        else
            error("GeoTIFF reading requires ArchGDAL.jl. Install with: using Pkg; Pkg.add(\"ArchGDAL\")")
        end
    elseif endswith(path, ".txt") || endswith(path, ".asc")
        # Simple text format fallback
        data = read_array(path)
        metadata = Dict{String,Any}(
            "crs" => "EPSG:4326",
            "dx" => 1.0,
            "dy" => 1.0,
            "x_origin" => 0.0,
            "y_origin" => 0.0,
        )
        return data, metadata
    else
        error("Unsupported file format: $path")
    end
end

"""
    read_scenario_toml(path::String)

Read a scenario configuration from TOML file.

# Expected TOML structure:
```toml
[scenario]
name = "My Scenario"
description = "Description"

[grid]
nx = 100
ny = 100
dx = 10.0

[input]
dem = "dem.txt"
roughness = 0.03
rainfall = "rainfall.csv"

[parameters]
t_end = 3600.0
dt_max = 1.0
cfl = 0.7

[output]
directory = "results"
interval = 60.0
points = [[50, 50], [25, 25]]
```

# Returns
- `Dict`: Parsed TOML configuration
"""
function read_scenario_toml(path::String)
    isfile(path) || throw(ArgumentError("Scenario file not found: $path"))
    TOML.parsefile(path)
end

"""
    load_scenario_from_toml(path::String)

Load a complete Scenario from a TOML configuration file.

# Arguments
- `path`: Path to scenario.toml file

# Returns
- `Scenario`: Complete scenario ready for simulation
"""
function load_scenario_from_toml(path::String)
    config = read_scenario_toml(path)
    base_dir = dirname(path)

    # Parse grid
    grid_config = config["grid"]
    nx = grid_config["nx"]
    ny = grid_config["ny"]
    dx = Float64(grid_config["dx"])
    dy = Float64(get(grid_config, "dy", dx))
    x_origin = Float64(get(grid_config, "x_origin", 0.0))
    y_origin = Float64(get(grid_config, "y_origin", 0.0))
    crs = get(grid_config, "crs", "EPSG:4326")

    grid = Grid(nx, ny, dx, dy, x_origin, y_origin, crs)

    # Parse input
    input_config = config["input"]

    # Load DEM
    dem_path = joinpath(base_dir, input_config["dem"])
    elevation = read_array(dem_path)
    size(elevation) == (nx, ny) || throw(DimensionMismatch(
        "DEM dimensions $(size(elevation)) don't match grid ($nx, $ny)"))

    # Load or create roughness
    roughness_input = input_config["roughness"]
    if roughness_input isa Number
        roughness = fill(Float64(roughness_input), nx, ny)
    else
        roughness_path = joinpath(base_dir, roughness_input)
        roughness = read_array(roughness_path)
    end

    topo = Topography(elevation, roughness, grid)

    # Load rainfall
    rainfall_path = joinpath(base_dir, input_config["rainfall"])
    rainfall = read_rainfall_csv(rainfall_path)

    # Parse parameters
    param_config = config["parameters"]
    params = SimulationParameters(
        dt_max = Float64(get(param_config, "dt_max", 1.0)),
        t_end = Float64(param_config["t_end"]),
        cfl = Float64(get(param_config, "cfl", 0.7)),
        h_min = Float64(get(param_config, "h_min", 0.001)),
        output_interval = Float64(get(param_config, "output_interval", 60.0)),
    )

    # Parse output
    output_config = get(config, "output", Dict())
    output_dir = joinpath(base_dir, get(output_config, "directory", "results"))
    points_raw = get(output_config, "points", Vector{Any}())
    output_points = Tuple{Int,Int}[(Int(p[1]), Int(p[2])) for p in points_raw]

    # Create scenario
    scenario_config = get(config, "scenario", Dict())
    name = get(scenario_config, "name", basename(dirname(path)))

    Scenario(name, grid, topo, params, rainfall, output_points, output_dir)
end

"""
    read_points_geojson(path::String)

Read output points from GeoJSON file.

# Returns
- `Vector{Tuple{Float64,Float64}}`: List of (x, y) coordinates
"""
function read_points_geojson(path::String)
    # Placeholder for GeoJSON support
    error("GeoJSON reading not yet implemented. Use TOML output points instead.")
end
