# Running Simulations

## Basic Usage

### From TOML Configuration

```julia
using HydroForge

# Run simulation from scenario file
results = HydroForge.run("path/to/scenario.toml")

# Access results
println("Run ID: ", results["run_id"])
println("Status: ", results["status"])
println("Max depth: ", results["max_depth"], " m")
println("Output: ", results["output_dir"])
```

### Demo Scenario

```julia
# Run bundled demo
HydroForge.run_demo()
```

## Programmatic Setup

For more control, build scenarios programmatically:

```julia
using HydroForge

# Create grid
grid = Grid(100, 100, 10.0)  # 100x100 cells, 10m resolution

# Create elevation (example: sloped plane)
elevation = zeros(100, 100)
for i in 1:100, j in 1:100
    elevation[i, j] = 10.0 + 0.01 * i + 0.01 * j
end

# Create topography
topo = Topography(elevation, 0.03, grid)  # Uniform Manning's n = 0.03

# Create rainfall event
times = [0.0, 1800.0, 3600.0, 5400.0, 7200.0]
intensities = [0.0, 30.0, 50.0, 20.0, 0.0]
rainfall = RainfallEvent(times, intensities)

# Create parameters
params = SimulationParameters(
    t_end = 7200.0,      # 2 hours
    dt_max = 1.0,
    cfl = 0.7,
    h_min = 0.001,
    output_interval = 300.0
)

# Create scenario
scenario = Scenario(
    "My Scenario",
    grid,
    topo,
    params,
    rainfall,
    [(50, 50), (25, 25)],  # Output points
    "results"
)

# Initialize state
state = SimulationState(grid)

# Run simulation
results = run_simulation!(state, scenario)
```

## Monitoring Progress

### Progress Callback

```julia
function my_progress(fraction)
    println("Progress: ", round(fraction * 100, digits=1), "%")
end

results = run_simulation!(state, scenario; progress_callback=my_progress)
```

### Accessing State During Simulation

```julia
# Manual stepping for custom control
workspace = SimulationWorkspace(grid)
state = SimulationState(grid)

while state.t < params.t_end
    dt = compute_dt(state, grid, params)
    step!(state, topo, params, rainfall, dt, workspace)

    # Access current state
    println("Time: ", state.t, " Max depth: ", max_depth(state))
end
```

## Performance Tips

### Memory

- Use `Float32` for large domains: `Grid{Float32}(...)`
- Preallocate workspace: `SimulationWorkspace(grid)`

### Speed

- Increase `dt_max` if stable (check CFL warnings)
- Reduce `output_interval` for fewer file writes
- Use smaller `h_min` only if needed for accuracy

### Typical Runtimes

| Grid Size | Simulation Time | Approx. Runtime |
|-----------|-----------------|-----------------|
| 64×64 | 1 hour | ~5 seconds |
| 256×256 | 1 hour | ~30 seconds |
| 1024×1024 | 1 hour | ~5 minutes |

*Benchmarks on Apple M1, single-threaded*

## Error Handling

```julia
try
    results = HydroForge.run("scenario.toml")
catch e
    if e isa ArgumentError
        println("Invalid input: ", e.msg)
    elseif e isa DimensionMismatch
        println("Data dimensions don't match: ", e.msg)
    else
        rethrow(e)
    end
end
```
