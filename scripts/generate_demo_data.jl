#!/usr/bin/env julia
# Generate demo data for HydroForge
# Creates synthetic DEM, roughness, and rainfall data

using DelimitedFiles

# Configuration
const NX = 64
const NY = 64
const DX = 10.0  # 10m resolution

# Output directory
const DEMO_DIR = joinpath(@__DIR__, "..", "assets", "demo_data")

"""
Generate a synthetic DEM with valley terrain.
- Higher elevation at edges
- Lower in center (drainage point)
- Some urban features (raised blocks)
"""
function generate_dem()
    elevation = zeros(NX, NY)

    # Base terrain: bowl shape (higher at edges)
    for i in 1:NX, j in 1:NY
        # Distance from center
        cx, cy = NX/2, NY/2
        dist = sqrt((i - cx)^2 + (j - cy)^2)

        # Bowl shape with gentle slope toward center
        elevation[i, j] = 10.0 + 0.1 * dist

        # Add some random micro-topography
        elevation[i, j] += 0.05 * randn()
    end

    # Add raised blocks to represent buildings
    # Create a grid of blocks
    block_size = 8
    block_height = 2.0

    for bi in 1:4, bj in 1:4
        # Skip some blocks for streets
        if (bi + bj) % 2 == 0
            continue
        end

        i_start = 8 + (bi - 1) * 12
        j_start = 8 + (bj - 1) * 12

        for i in i_start:min(i_start + block_size - 1, NX)
            for j in j_start:min(j_start + block_size - 1, NY)
                elevation[i, j] += block_height
            end
        end
    end

    # Create a drainage channel toward outlet (lower-right)
    for i in Int(NX/2):NX
        j = Int(NY/2)
        for dj in -2:2
            if 1 <= j + dj <= NY
                elevation[i, j + dj] -= 0.5
            end
        end
    end

    return elevation
end

"""
Generate roughness map with urban pattern.
- Streets: n = 0.013 (smooth concrete)
- Vegetated areas: n = 0.035 (grass)
- Buildings: n = 0.10 (high roughness for obstacles)
"""
function generate_roughness(elevation)
    roughness = fill(0.03, NX, NY)  # Default grass

    # Streets (lower roughness) - grid pattern
    for i in 1:NX
        for j in 1:NY
            # Vertical streets every 12 cells
            if (i - 4) % 12 < 4
                roughness[i, j] = 0.013
            end
            # Horizontal streets every 12 cells
            if (j - 4) % 12 < 4
                roughness[i, j] = 0.013
            end
        end
    end

    # Buildings (high roughness) - where elevation is raised
    base_elev = 10.0
    for i in 1:NX, j in 1:NY
        if elevation[i, j] > base_elev + 1.5
            roughness[i, j] = 0.10
        end
    end

    return roughness
end

"""
Generate rainfall CSV files.
"""
function generate_rainfall()
    # Short intense storm (1 hour, 50 mm/hr peak)
    short_times = [0.0, 600.0, 1200.0, 1800.0, 2400.0, 3000.0, 3600.0]
    short_intensities = [0.0, 20.0, 50.0, 50.0, 30.0, 10.0, 0.0]

    # Long moderate storm (3 hours)
    long_times = [0.0, 1800.0, 3600.0, 5400.0, 7200.0, 9000.0, 10800.0]
    long_intensities = [0.0, 15.0, 25.0, 20.0, 15.0, 5.0, 0.0]

    return (
        short = (times=short_times, intensities=short_intensities),
        long = (times=long_times, intensities=long_intensities)
    )
end

"""
Write all demo data files.
"""
function generate_all()
    mkpath(DEMO_DIR)

    println("Generating demo data in $DEMO_DIR")

    # Generate and write DEM
    println("  - Generating DEM...")
    dem = generate_dem()
    writedlm(joinpath(DEMO_DIR, "dem.txt"), dem, ' ')

    # Generate and write roughness
    println("  - Generating roughness map...")
    roughness = generate_roughness(dem)
    writedlm(joinpath(DEMO_DIR, "roughness.txt"), roughness, ' ')

    # Generate and write rainfall
    println("  - Generating rainfall data...")
    rainfall = generate_rainfall()

    open(joinpath(DEMO_DIR, "rainfall_short.csv"), "w") do io
        println(io, "time,intensity")
        for (t, i) in zip(rainfall.short.times, rainfall.short.intensities)
            println(io, "$t,$i")
        end
    end

    open(joinpath(DEMO_DIR, "rainfall_long.csv"), "w") do io
        println(io, "time,intensity")
        for (t, i) in zip(rainfall.long.times, rainfall.long.intensities)
            println(io, "$t,$i")
        end
    end

    # Write scenario TOML
    println("  - Writing scenario configuration...")
    scenario_toml = """
[scenario]
name = "Minimal City Demo"
description = "64x64 synthetic urban catchment for demonstration"

[grid]
nx = $NX
ny = $NY
dx = $DX
x_origin = 0.0
y_origin = 0.0
crs = "EPSG:32632"

[input]
dem = "dem.txt"
roughness = "roughness.txt"
rainfall = "rainfall_short.csv"

[parameters]
t_end = 3600.0
dt_max = 1.0
cfl = 0.7
h_min = 0.001
output_interval = 60.0

[output]
directory = "results"
points = [[32, 32], [16, 16], [48, 48]]
"""

    open(joinpath(DEMO_DIR, "scenario.toml"), "w") do io
        write(io, scenario_toml)
    end

    println("Demo data generation complete!")
    println("Files created:")
    for f in readdir(DEMO_DIR)
        println("  - $f")
    end
end

# Run if executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    generate_all()
end
