# HydroForge Benchmarks
# Performance benchmarks for core solver components

using BenchmarkTools
using HydroForge

const SUITE = BenchmarkGroup()

# =============================================================================
# Grid Operations
# =============================================================================
SUITE["grid"] = BenchmarkGroup()

SUITE["grid"]["construction_small"] = @benchmarkable Grid(64, 64, 10.0)
SUITE["grid"]["construction_medium"] = @benchmarkable Grid(256, 256, 10.0)
SUITE["grid"]["construction_large"] = @benchmarkable Grid(1024, 1024, 10.0)

# =============================================================================
# State Operations
# =============================================================================
SUITE["state"] = BenchmarkGroup()

function setup_state(n)
    grid = Grid(n, n, 10.0)
    state = SimulationState(grid)
    # Add some water
    for i in 1:n, j in 1:n
        state.h[i, j] = rand() * 0.1
    end
    return state, grid
end

for n in [64, 256, 1024]
    SUITE["state"]["total_volume_$(n)"] = @benchmarkable total_volume(s, g) setup=(s, g = setup_state($n))
    SUITE["state"]["max_depth_$(n)"] = @benchmarkable max_depth(s) setup=(s, _ = setup_state($n))
    SUITE["state"]["wet_cells_$(n)"] = @benchmarkable wet_cells(s, 0.001) setup=(s, _ = setup_state($n))
end

# =============================================================================
# Timestep Computation
# =============================================================================
SUITE["timestep"] = BenchmarkGroup()

function setup_timestep(n)
    grid = Grid(n, n, 10.0)
    state = SimulationState(grid)
    params = SimulationParameters()
    # Add water with varying depths
    for i in 1:n, j in 1:n
        state.h[i, j] = rand() * 0.5
        state.qx[i, j] = rand() * 0.01
        state.qy[i, j] = rand() * 0.01
    end
    return state, grid, params
end

for n in [64, 256, 1024]
    SUITE["timestep"]["compute_dt_$(n)"] = @benchmarkable compute_dt(s, g, p) setup=(s, g, p = setup_timestep($n))
end

# =============================================================================
# Flux Computation
# =============================================================================
SUITE["flux"] = BenchmarkGroup()

function setup_flux(n)
    grid = Grid(n, n, 10.0)
    state = SimulationState(grid)
    elevation = zeros(n, n)
    for i in 1:n, j in 1:n
        elevation[i, j] = (i + j) * 0.01  # Gentle slope
        state.h[i, j] = rand() * 0.2
    end
    topo = Topography(elevation, 0.03, grid)
    params = SimulationParameters()
    return state, topo, grid, params
end

for n in [64, 256]
    SUITE["flux"]["compute_flux_x_$(n)"] = @benchmarkable begin
        compute_flux_x!(s.qx, s, t, g, p, 0.1)
    end setup=(s, t, g, p = setup_flux($n))

    SUITE["flux"]["compute_flux_y_$(n)"] = @benchmarkable begin
        compute_flux_y!(s.qy, s, t, g, p, 0.1)
    end setup=(s, t, g, p = setup_flux($n))
end

# =============================================================================
# Full Timestep
# =============================================================================
SUITE["solver"] = BenchmarkGroup()

function setup_solver(n)
    grid = Grid(n, n, 10.0)
    state = SimulationState(grid)
    elevation = zeros(n, n)
    for i in 1:n, j in 1:n
        elevation[i, j] = (i + j) * 0.01
        state.h[i, j] = rand() * 0.1
    end
    topo = Topography(elevation, 0.03, grid)
    params = SimulationParameters(dt_max=1.0)
    rainfall = RainfallEvent([0.0, 3600.0], [20.0, 0.0])

    workspace = SimulationWorkspace(grid)
    return state, topo, params, rainfall, workspace
end

for n in [64, 256]
    SUITE["solver"]["step_$(n)"] = @benchmarkable begin
        step!(ws, s, t, p, r, 0.1)
    end setup=(s, t, p, r, ws = setup_solver($n))
end

# =============================================================================
# Run Benchmarks
# =============================================================================
function run_benchmarks(; verbose=true)
    results = run(SUITE, verbose=verbose)
    return results
end

# Entry point for running benchmarks
if abspath(PROGRAM_FILE) == @__FILE__
    results = run_benchmarks()
    display(results)
end
