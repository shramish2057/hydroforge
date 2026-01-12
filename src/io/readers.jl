# HydroForge IO Readers
# Functions for reading input data files

using DelimitedFiles
using TOML
using JSON3

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

Supports multiple raster formats:
- `.tif`, `.tiff`: GeoTIFF (requires ArchGDAL.jl)
- `.asc`: ESRI ASCII Grid format
- `.txt`: Simple space-separated text format

# Returns
- `data::Matrix{Float64}`: Raster data
- `metadata::Dict`: Contains crs, dx, dy, x_origin, y_origin

# Note
For GeoTIFF support, install ArchGDAL.jl:
```julia
using Pkg; Pkg.add("ArchGDAL")
using ArchGDAL
```
"""
function read_geotiff(path::String)
    isfile(path) || throw(ArgumentError("Raster file not found: $path"))

    ext = lowercase(splitext(path)[2])

    if ext in [".tif", ".tiff"]
        # Try to use ArchGDAL if available
        if isdefined(Main, :ArchGDAL)
            return _read_geotiff_archgdal(path)
        else
            error("""GeoTIFF reading requires ArchGDAL.jl.

Install with:
    using Pkg; Pkg.add("ArchGDAL")

Then load before HydroForge:
    using ArchGDAL
    using HydroForge

Alternatively, convert your GeoTIFF to ASCII format (.asc) which is natively supported.""")
        end
    elseif ext == ".asc"
        # ESRI ASCII Grid format
        return read_esri_ascii(path)
    elseif ext == ".txt"
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
        error("Unsupported raster format: $ext. Supported: .tif, .tiff, .asc, .txt")
    end
end

"""
    _read_geotiff_archgdal(path::String)

Internal function to read GeoTIFF using ArchGDAL.
Only called if ArchGDAL is available.
"""
function _read_geotiff_archgdal(path::String)
    # This uses Main.ArchGDAL to avoid compile-time dependency
    AG = Main.ArchGDAL

    dataset = AG.read(path)
    band = AG.getband(dataset, 1)

    # Read data
    data = Matrix{Float64}(AG.read(band))

    # Get geotransform: [x_origin, dx, rot1, y_origin, rot2, -dy]
    gt = AG.getgeotransform(dataset)

    # Get projection
    crs = try
        proj = AG.getproj(dataset)
        if isempty(proj)
            "EPSG:4326"
        else
            proj
        end
    catch
        "EPSG:4326"
    end

    metadata = Dict{String,Any}(
        "crs" => crs,
        "dx" => abs(gt[2]),
        "dy" => abs(gt[6]),
        "x_origin" => gt[1],
        "y_origin" => gt[4] + gt[6] * size(data, 1),  # Adjust for top-left origin
        "nodata" => try AG.getnodatavalue(band) catch; nothing end,
    )

    # Handle nodata values
    nodata = metadata["nodata"]
    if nodata !== nothing
        data[data .== nodata] .= NaN
    end

    return data, metadata
end

"""
    read_esri_ascii(path::String)

Read an ESRI ASCII Grid file (.asc format).

# Expected format:
```
ncols         100
nrows         100
xllcorner     0.0
yllcorner     0.0
cellsize      10.0
NODATA_value  -9999
1.0 2.0 3.0 ...
```

# Returns
- `data::Matrix{Float64}`: Raster data
- `metadata::Dict`: Contains crs, dx, dy, x_origin, y_origin, nodata
"""
function read_esri_ascii(path::String)
    isfile(path) || throw(ArgumentError("ASCII grid file not found: $path"))

    # Read file
    lines = readlines(path)

    # Parse header
    header = Dict{String,Any}()
    data_start = 1

    for (i, line) in enumerate(lines)
        parts = split(strip(line))
        if length(parts) >= 2
            key = lowercase(parts[1])
            if key in ["ncols", "nrows"]
                header[key] = parse(Int, parts[2])
                data_start = i + 1
            elseif key in ["xllcorner", "xllcenter", "yllcorner", "yllcenter",
                          "cellsize", "dx", "dy", "nodata_value"]
                header[key] = parse(Float64, parts[2])
                data_start = i + 1
            else
                # Not a header line, data starts here
                break
            end
        else
            break
        end
    end

    # Extract dimensions
    ncols = get(header, "ncols", nothing)
    nrows = get(header, "nrows", nothing)

    if ncols === nothing || nrows === nothing
        error("ASCII grid missing ncols or nrows header")
    end

    # Extract cell size
    cellsize = get(header, "cellsize", 1.0)
    dx = get(header, "dx", cellsize)
    dy = get(header, "dy", cellsize)

    # Extract origin (handle both corner and center)
    if haskey(header, "xllcorner")
        x_origin = header["xllcorner"]
    elseif haskey(header, "xllcenter")
        x_origin = header["xllcenter"] - dx / 2
    else
        x_origin = 0.0
    end

    if haskey(header, "yllcorner")
        y_origin = header["yllcorner"]
    elseif haskey(header, "yllcenter")
        y_origin = header["yllcenter"] - dy / 2
    else
        y_origin = 0.0
    end

    nodata = get(header, "nodata_value", -9999.0)

    # Read data
    data = Matrix{Float64}(undef, nrows, ncols)
    row = 1

    for i in data_start:length(lines)
        if row > nrows
            break
        end

        values = parse.(Float64, split(strip(lines[i])))
        if length(values) == ncols
            data[row, :] = values
            row += 1
        elseif length(values) > 0
            # Handle multiple values per line that might span multiple rows
            for v in values
                if row <= nrows
                    col = (row - 1) * ncols + 1
                    # ... simplified: assume one row per line
                end
            end
        end
    end

    # Handle nodata
    if nodata !== nothing
        data[data .== nodata] .= NaN
    end

    metadata = Dict{String,Any}(
        "crs" => "EPSG:4326",  # ASCII grids don't store CRS
        "dx" => dx,
        "dy" => dy,
        "x_origin" => x_origin,
        "y_origin" => y_origin,
        "nodata" => nodata,
        "ncols" => ncols,
        "nrows" => nrows,
    )

    return data, metadata
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

Supports both Point and MultiPoint geometries from a FeatureCollection.

# Expected GeoJSON format:
```json
{
  "type": "FeatureCollection",
  "features": [
    {
      "type": "Feature",
      "geometry": {
        "type": "Point",
        "coordinates": [100.5, 50.5]
      },
      "properties": {"name": "Gauge 1"}
    }
  ]
}
```

# Returns
- `Vector{Tuple{Float64,Float64}}`: List of (x, y) coordinates

# Note
Coordinates are returned as (x, y) = (longitude, latitude) or (easting, northing).
For grid cell indices, use `point_to_cell()` to convert.
"""
function read_points_geojson(path::String)
    isfile(path) || throw(ArgumentError("GeoJSON file not found: $path"))

    # Read and parse JSON
    json_str = read(path, String)
    geojson = JSON3.read(json_str)

    points = Tuple{Float64,Float64}[]

    # Handle different GeoJSON structures
    if haskey(geojson, :type)
        geotype = geojson[:type]

        if geotype == "FeatureCollection"
            # FeatureCollection with multiple features
            for feature in geojson[:features]
                _extract_points!(points, feature[:geometry])
            end
        elseif geotype == "Feature"
            # Single feature
            _extract_points!(points, geojson[:geometry])
        elseif geotype == "GeometryCollection"
            # GeometryCollection
            for geom in geojson[:geometries]
                _extract_points!(points, geom)
            end
        else
            # Direct geometry
            _extract_points!(points, geojson)
        end
    else
        throw(ArgumentError("Invalid GeoJSON: missing 'type' field"))
    end

    if isempty(points)
        @warn "No points found in GeoJSON file: $path"
    end

    points
end

"""
    _extract_points!(points, geometry)

Extract point coordinates from a GeoJSON geometry object (internal helper).
"""
function _extract_points!(points::Vector{Tuple{Float64,Float64}}, geometry)
    geomtype = geometry[:type]

    if geomtype == "Point"
        coords = geometry[:coordinates]
        push!(points, (Float64(coords[1]), Float64(coords[2])))
    elseif geomtype == "MultiPoint"
        for coord in geometry[:coordinates]
            push!(points, (Float64(coord[1]), Float64(coord[2])))
        end
    elseif geomtype == "LineString"
        # Extract all vertices as points
        for coord in geometry[:coordinates]
            push!(points, (Float64(coord[1]), Float64(coord[2])))
        end
    elseif geomtype == "Polygon"
        # Extract first ring's vertices (outer boundary)
        for coord in geometry[:coordinates][1]
            push!(points, (Float64(coord[1]), Float64(coord[2])))
        end
    end
    # Silently ignore other geometry types (MultiPolygon, etc.)

    nothing
end

"""
    point_to_cell(x::Real, y::Real, grid::Grid)

Convert world coordinates (x, y) to grid cell indices (i, j).

# Returns
- `Tuple{Int,Int}`: Cell indices (i, j) where 1 ≤ i ≤ nx, 1 ≤ j ≤ ny

# Note
Returns indices clamped to grid bounds if point is outside grid.
"""
function point_to_cell(x::Real, y::Real, grid::Grid{T}) where T
    i = Int(floor((x - grid.x_origin) / grid.dx)) + 1
    j = Int(floor((y - grid.y_origin) / grid.dy)) + 1

    # Clamp to grid bounds
    i = clamp(i, 1, grid.nx)
    j = clamp(j, 1, grid.ny)

    (i, j)
end

"""
    read_points_geojson(path::String, grid::Grid)

Read output points from GeoJSON and convert to grid cell indices.

# Returns
- `Vector{Tuple{Int,Int}}`: List of (i, j) cell indices
"""
function read_points_geojson(path::String, grid::Grid)
    coords = read_points_geojson(path)
    [point_to_cell(x, y, grid) for (x, y) in coords]
end

"""
    read_network_toml(path::String)

Read a drainage network configuration from TOML file.

# Expected TOML structure:
```toml
[network]
name = "My Network"
description = "Description"

[[junctions]]
id = 1
x = 100.0
y = 100.0
invert = -2.0
ground = 0.0
type = "MANHOLE"  # MANHOLE, STORAGE, OUTFALL, DIVIDER
init_depth = 0.0
ponded_area = 0.0

[[pipes]]
id = 1
upstream_node = 1
downstream_node = 2
shape = "circular"  # circular or rectangular
diameter = 0.45     # for circular pipes
# width = 1.0       # for rectangular pipes
# height = 0.5      # for rectangular pipes
length = 100.0
roughness = 0.013
invert_up = -2.0
invert_down = -2.5

[[inlets]]
id = 1
junction_id = 1
grid_i = 10
grid_j = 10
type = "GRATE"  # GRATE, CURB, COMBINATION, SLOTTED, DROP
length = 0.6
width = 0.6
clogging_factor = 1.0

[[outlets]]
id = 1
junction_id = 2
type = "FREE"  # FREE, FIXED, TIDAL, FLAP
invert = -2.5
fixed_stage = 0.0
```

# Returns
- `Dict`: Parsed TOML configuration
"""
function read_network_toml(path::String)
    isfile(path) || throw(ArgumentError("Network file not found: $path"))
    TOML.parsefile(path)
end

"""
    load_network_from_toml(path::String)

Load a complete DrainageNetwork from a TOML configuration file.

# Arguments
- `path`: Path to network.toml file

# Returns
- `DrainageNetwork`: Complete drainage network ready for simulation
"""
function load_network_from_toml(path::String)
    config = read_network_toml(path)

    # Parse junctions
    junctions = Junction{Float64}[]
    junction_configs = get(config, "junctions", [])

    for jc in junction_configs
        id = jc["id"]
        x = Float64(jc["x"])
        y = Float64(jc["y"])
        invert = Float64(jc["invert"])
        ground = Float64(jc["ground"])

        # Parse junction type
        type_str = uppercase(get(jc, "type", "MANHOLE"))
        junction_type = if type_str == "MANHOLE"
            MANHOLE
        elseif type_str == "STORAGE"
            STORAGE
        elseif type_str == "OUTFALL"
            OUTFALL
        elseif type_str == "DIVIDER"
            DIVIDER
        else
            @warn "Unknown junction type '$type_str', defaulting to MANHOLE"
            MANHOLE
        end

        # Optional fields
        init_depth = Float64(get(jc, "init_depth", 0.0))
        ponded_area = Float64(get(jc, "ponded_area", 0.0))
        max_depth = Float64(get(jc, "max_depth", ground - invert))

        # Storage curve (optional)
        storage_curve = nothing
        if haskey(jc, "storage_curve")
            curve_data = jc["storage_curve"]
            storage_curve = Tuple{Float64,Float64}[(Float64(p[1]), Float64(p[2])) for p in curve_data]
        end

        junction = Junction{Float64}(
            id, x, y, invert, ground, junction_type,
            storage_curve, max_depth, init_depth, ponded_area
        )
        push!(junctions, junction)
    end

    # Parse pipes
    pipes = PipeSegment{Float64}[]
    pipe_configs = get(config, "pipes", [])

    for pc in pipe_configs
        id = pc["id"]
        upstream_node = pc["upstream_node"]
        downstream_node = pc["downstream_node"]
        length = Float64(pc["length"])
        roughness = Float64(get(pc, "roughness", 0.013))
        invert_up = Float64(pc["invert_up"])
        invert_down = Float64(pc["invert_down"])

        # Parse cross-section
        shape = lowercase(get(pc, "shape", "circular"))
        section::PipeCrossSection = if shape == "circular"
            diameter = Float64(pc["diameter"])
            CircularPipe(diameter)
        elseif shape == "rectangular"
            width = Float64(pc["width"])
            height = Float64(pc["height"])
            RectangularPipe(width, height)
        else
            error("Unknown pipe shape: $shape. Use 'circular' or 'rectangular'.")
        end

        pipe = PipeSegment{Float64}(
            id, upstream_node, downstream_node, section,
            length, roughness, invert_up, invert_down
        )
        push!(pipes, pipe)
    end

    # Parse inlets
    inlets = Inlet{Float64}[]
    inlet_configs = get(config, "inlets", [])

    for ic in inlet_configs
        id = ic["id"]
        junction_id = ic["junction_id"]
        grid_i = ic["grid_i"]
        grid_j = ic["grid_j"]

        # Parse inlet type
        type_str = uppercase(get(ic, "type", "GRATE"))
        inlet_type = if type_str == "GRATE"
            GRATE
        elseif type_str == "CURB"
            CURB
        elseif type_str == "COMBINATION"
            COMBINATION
        elseif type_str == "SLOTTED"
            SLOTTED
        elseif type_str == "DROP"
            DROP
        else
            @warn "Unknown inlet type '$type_str', defaulting to GRATE"
            GRATE
        end

        # Optional fields with defaults
        length = Float64(get(ic, "length", 0.6))
        width = Float64(get(ic, "width", 0.6))
        clogging_factor = Float64(get(ic, "clogging_factor", 1.0))
        weir_coeff = Float64(get(ic, "weir_coeff", 1.66))
        orifice_coeff = Float64(get(ic, "orifice_coeff", 0.67))
        depression = Float64(get(ic, "depression", 0.0))

        inlet = Inlet{Float64}(
            id, junction_id, grid_i, grid_j, inlet_type,
            length, width, clogging_factor, weir_coeff, orifice_coeff, depression
        )
        push!(inlets, inlet)
    end

    # Parse outlets
    outlets = Outlet{Float64}[]
    outlet_configs = get(config, "outlets", [])

    for oc in outlet_configs
        id = oc["id"]
        junction_id = oc["junction_id"]

        # Optional fields
        grid_i = get(oc, "grid_i", -1)
        grid_j = get(oc, "grid_j", -1)
        invert = Float64(get(oc, "invert", 0.0))
        fixed_stage = Float64(get(oc, "fixed_stage", 0.0))
        flap_loss = Float64(get(oc, "flap_loss", 0.0))

        # Parse outlet type
        type_str = uppercase(get(oc, "type", "FREE"))
        outlet_type = Symbol(type_str)

        # Tide curve (optional)
        tide_curve = nothing
        if haskey(oc, "tide_curve")
            curve_data = oc["tide_curve"]
            tide_curve = Tuple{Float64,Float64}[(Float64(p[1]), Float64(p[2])) for p in curve_data]
        end

        outlet = Outlet{Float64}(
            id, junction_id, grid_i, grid_j, outlet_type,
            invert, fixed_stage, tide_curve, flap_loss
        )
        push!(outlets, outlet)
    end

    # Validate we have at least one junction
    if isempty(junctions)
        throw(ArgumentError("Network must have at least one junction"))
    end

    # Create network
    network = DrainageNetwork(pipes, junctions, inlets, outlets)

    # Validate network
    issues = validate(network)
    if !isempty(issues)
        @warn "Network validation found $(length(issues)) issue(s):" issues
    end

    network
end
