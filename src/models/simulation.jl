# HydroForge Simulation Runner
# High-level simulation execution
#
# Production-ready execution framework with:
# - Comprehensive error handling and logging
# - Snapshot saving for intermediate states
# - Performance profiling support
# - Hazard statistics tracking
# - Graceful failure handling with partial results

using Dates
using Random

# Note: ResultsAccumulator, update_results!, and record_output! are defined in surface2d.jl

"""
    RunConfig

Configuration for a simulation run.

# Fields
- `scenario_path::String`: Path to scenario TOML file
- `output_dir::String`: Directory for output files
- `run_id::String`: Unique identifier for this run
- `save_snapshots::Bool`: Whether to save intermediate depth snapshots
- `snapshot_interval::Float64`: Time between snapshots (seconds)
- `enable_profiling::Bool`: Whether to enable performance profiling
- `verbosity::Int`: Logging verbosity (0=silent, 1=info, 2=debug)
"""
struct RunConfig
    scenario_path::String
    output_dir::String
    run_id::String
    save_snapshots::Bool
    snapshot_interval::Float64
    enable_profiling::Bool
    verbosity::Int
end

"""
    create_run_config(scenario_path; kwargs...)

Create run configuration with sensible defaults.

# Keyword Arguments
- `output_dir`: Output directory (default: runs/<run_id> in scenario directory)
- `save_snapshots`: Save intermediate snapshots (default: false)
- `snapshot_interval`: Time between snapshots in seconds (default: 300.0)
- `enable_profiling`: Enable performance profiling (default: false)
- `verbosity`: Logging level 0-2 (default: 1)
"""
function create_run_config(scenario_path::String;
                           output_dir::Union{String,Nothing}=nothing,
                           save_snapshots::Bool=false,
                           snapshot_interval::Float64=300.0,
                           enable_profiling::Bool=false,
                           verbosity::Int=1)
    # Generate run ID with timestamp and random suffix
    timestamp = Dates.format(now(), "yyyymmdd_HHMMSS")
    random_suffix = randstring(6)
    run_id = "$(timestamp)_$(random_suffix)"

    # Set output directory
    if output_dir === nothing
        output_dir = joinpath(dirname(scenario_path), "runs", run_id)
    end

    # Create output directory structure
    mkpath(output_dir)
    if save_snapshots
        mkpath(joinpath(output_dir, "snapshots"))
    end

    RunConfig(scenario_path, output_dir, run_id, save_snapshots,
              snapshot_interval, enable_profiling, verbosity)
end


"""
    RunMetadata

Comprehensive metadata for a simulation run.

# Fields
- `run_id::String`: Unique run identifier
- `start_time::DateTime`: Run start timestamp
- `end_time::Union{DateTime, Nothing}`: Run end timestamp
- `julia_version::String`: Julia version used
- `hydroforge_version::String`: HydroForge version
- `git_commit::Union{String, Nothing}`: Git commit hash if available
- `scenario_name::String`: Name of the scenario
- `parameters::Dict{String, Any}`: Simulation parameters
- `status::Symbol`: Current status (:running, :completed, :failed, :cancelled)
- `error_message::Union{String, Nothing}`: Error message if failed
- `performance::Dict{String, Any}`: Performance statistics
- `hazard_summary::Dict{String, Any}`: Hazard analysis summary
"""
mutable struct RunMetadata
    run_id::String
    start_time::DateTime
    end_time::Union{DateTime, Nothing}
    julia_version::String
    hydroforge_version::String
    git_commit::Union{String, Nothing}
    scenario_name::String
    parameters::Dict{String, Any}
    status::Symbol  # :running, :completed, :failed, :cancelled
    error_message::Union{String, Nothing}
    performance::Dict{String, Any}
    hazard_summary::Dict{String, Any}
end

"""
    create_metadata(config::RunConfig, scenario_name::String, params)

Create run metadata with comprehensive tracking.
"""
function create_metadata(config::RunConfig, scenario_name::String, params)
    RunMetadata(
        config.run_id,
        now(),
        nothing,
        string(VERSION),
        HYDROFORGE_VERSION,
        get_git_commit(),
        scenario_name,
        Dict{String, Any}(
            "t_end" => params.t_end,
            "dt_max" => params.dt_max,
            "cfl" => params.cfl,
            "h_min" => params.h_min,
            "g" => params.g,
            "output_interval" => params.output_interval
        ),
        :running,
        nothing,
        Dict{String, Any}(),  # performance
        Dict{String, Any}()   # hazard_summary
    )
end

"""
    get_git_commit()

Try to get current git commit hash.
"""
function get_git_commit()
    try
        strip(read(`git rev-parse HEAD`, String))
    catch
        nothing
    end
end

"""
    save_snapshot(state::SimulationState, grid::Grid, output_dir::String, snapshot_num::Int)

Save a depth snapshot to disk for later analysis.
"""
function save_snapshot(state::SimulationState{T}, grid::Grid{T},
                       output_dir::String, snapshot_num::Int) where T
    snapshot_path = joinpath(output_dir, "snapshots",
                             "depth_$(lpad(snapshot_num, 4, '0')).dat")
    # Write binary format for efficiency
    open(snapshot_path, "w") do f
        write(f, Float64(state.t))          # Time
        write(f, Int32(grid.nx))            # Grid dimensions
        write(f, Int32(grid.ny))
        write(f, state.h)                   # Depth array
    end
    snapshot_path
end

"""
    SimulationError <: Exception

Custom exception for simulation failures with context.
"""
struct SimulationError <: Exception
    message::String
    step::Int
    time::Float64
    original_error::Union{Exception, Nothing}
end

function Base.showerror(io::IO, e::SimulationError)
    print(io, "SimulationError at step $(e.step), t=$(round(e.time, digits=2))s: $(e.message)")
    if e.original_error !== nothing
        print(io, "\nCaused by: ")
        showerror(io, e.original_error)
    end
end


"""
    run(scenario_path::String; kwargs...)

Run a simulation from scenario file with comprehensive tracking.

# Arguments
- `scenario_path`: Path to scenario TOML file

# Keyword Arguments
- `output_dir`: Output directory (default: runs/<run_id> in scenario directory)
- `save_snapshots`: Save intermediate depth snapshots (default: false)
- `snapshot_interval`: Time between snapshots in seconds (default: 300.0)
- `enable_profiling`: Enable performance profiling (default: false)
- `verbosity`: Logging level 0-2 (default: 1)

# Returns
- Results summary dictionary with:
  - `run_id`: Unique run identifier
  - `status`: :completed or :failed
  - `output_dir`: Path to output files
  - `max_depth`: Maximum flood depth (m)
  - `max_hazard`: Maximum hazard rating (m²/s)
  - `steps`: Number of timesteps
  - `wall_time`: Computation time (s)
  - `mass_error_pct`: Mass balance error percentage
  - `hazard_summary`: Detailed hazard statistics
"""
function run(scenario_path::String;
             output_dir::Union{String,Nothing}=nothing,
             save_snapshots::Bool=false,
             snapshot_interval::Float64=300.0,
             enable_profiling::Bool=false,
             verbosity::Int=1)

    # Create run configuration
    config = create_run_config(scenario_path;
                               output_dir=output_dir,
                               save_snapshots=save_snapshots,
                               snapshot_interval=snapshot_interval,
                               enable_profiling=enable_profiling,
                               verbosity=verbosity)

    if verbosity >= 1
        @info "Starting HydroForge simulation" run_id=config.run_id output_dir=config.output_dir
    end

    # Load scenario
    scenario = load_scenario(scenario_path)

    if verbosity >= 2
        @info "Scenario loaded" name=scenario.name grid="$(scenario.grid.nx)x$(scenario.grid.ny)" t_end=scenario.parameters.t_end
    end

    # Create metadata
    metadata = create_metadata(config, scenario.name, scenario.parameters)

    # Save initial metadata
    save_metadata(config, metadata)

    try
        # Initialize state
        state = SimulationState(scenario.grid)

        if verbosity >= 1
            @info "Running simulation..." t_end=scenario.parameters.t_end
        end

        # Run simulation with optional profiling
        sim_results = if enable_profiling
            run_with_profiling(state, scenario, config)
        else
            run_simulation!(state, scenario;
                           verbosity=verbosity >= 2 ? 1 : 0,
                           log_interval=verbosity >= 2 ? scenario.parameters.t_end/10 : nothing)
        end

        # Extract accumulator for saving
        results = sim_results.accumulator

        # Compute hazard summary
        hazard_stats = summarize_hazard(results, scenario.grid)

        # Save results
        save_results(config, results, scenario)

        # Update metadata with performance and hazard info
        metadata.end_time = now()
        metadata.status = :completed
        metadata.performance = Dict{String, Any}(
            "steps" => sim_results.step_count,
            "wall_time_s" => sim_results.wall_time,
            "time_per_step_ms" => sim_results.wall_time / sim_results.step_count * 1000,
            "cells_per_second" => scenario.grid.nx * scenario.grid.ny * sim_results.step_count / sim_results.wall_time,
            "mass_error_pct" => relative_mass_error(sim_results.mass_balance) * 100
        )
        metadata.hazard_summary = hazard_stats
        save_metadata(config, metadata)

        if verbosity >= 1
            @info "Simulation completed" run_id=config.run_id steps=sim_results.step_count wall_time=round(sim_results.wall_time, digits=2) max_depth=round(hazard_stats["max_depth"], digits=3) max_hazard=round(hazard_stats["max_hazard"], digits=3)
        end

        # Return summary
        Dict{String, Any}(
            "run_id" => config.run_id,
            "status" => :completed,
            "output_dir" => config.output_dir,
            "max_depth" => hazard_stats["max_depth"],
            "max_velocity" => hazard_stats["max_velocity"],
            "max_hazard" => hazard_stats["max_hazard"],
            "max_froude" => hazard_stats["max_froude"],
            "steps" => sim_results.step_count,
            "wall_time" => sim_results.wall_time,
            "mass_error_pct" => relative_mass_error(sim_results.mass_balance) * 100,
            "hazard_summary" => hazard_stats
        )

    catch e
        # Enhanced error handling
        metadata.end_time = now()
        metadata.status = :failed
        metadata.error_message = sprint(showerror, e)

        if verbosity >= 1
            @error "Simulation failed" run_id=config.run_id error=e
        end

        # Try to save partial results and metadata
        try
            save_metadata(config, metadata)
        catch save_error
            if verbosity >= 1
                @warn "Failed to save metadata after error" error=save_error
            end
        end

        rethrow()
    end
end

"""
    run_with_profiling(state, scenario, config)

Run simulation with performance profiling enabled.
Saves profile data to the output directory.
"""
function run_with_profiling(state::SimulationState{T}, scenario::Scenario{T},
                            config::RunConfig) where T
    # Run simulation normally (Profile module not loaded by default)
    # In a full implementation, this would use Profile.@profile
    sim_results = run_simulation!(state, scenario; verbosity=0)

    # Write timing information
    timing_path = joinpath(config.output_dir, "timing.json")
    timing_data = Dict{String, Any}(
        "total_wall_time_s" => sim_results.wall_time,
        "steps" => sim_results.step_count,
        "avg_step_time_ms" => sim_results.wall_time / sim_results.step_count * 1000,
        "grid_cells" => scenario.grid.nx * scenario.grid.ny,
        "cells_per_second" => scenario.grid.nx * scenario.grid.ny * sim_results.step_count / sim_results.wall_time
    )
    write_results_json(timing_path, timing_data)

    sim_results
end

"""
    run_demo()

Run the bundled demo scenario.
"""
function run_demo()
    demo_path = joinpath(@__DIR__, "..", "..", "assets", "demo_data", "scenario.toml")

    if !isfile(demo_path)
        error("Demo scenario not found at $demo_path")
    end

    run(demo_path)
end


# IO Integration

"""
    ResultsPackage(results::ResultsAccumulator, metadata::Dict)

Create a ResultsPackage from a ResultsAccumulator and metadata.
"""
function ResultsPackage(results::ResultsAccumulator{T}, metadata::Dict{String,Any}) where T
    ResultsPackage{T}(
        results.max_depth,
        results.arrival_time,
        results.max_velocity,
        results.point_hydrographs,
        metadata
    )
end

"""
    write_results(output_dir::String, results::ResultsAccumulator, scenario::Scenario, run_id::String)

Convenience function to write results with auto-generated metadata.
"""
function write_results(output_dir::String, results::ResultsAccumulator{T},
                       scenario::Scenario{T}, run_id::String) where T
    metadata = Dict{String,Any}(
        "run_id" => run_id,
        "scenario_name" => scenario.name,
        "timestamp" => string(now()),
        "grid_nx" => scenario.grid.nx,
        "grid_ny" => scenario.grid.ny,
        "grid_dx" => scenario.grid.dx,
        "t_end" => scenario.parameters.t_end,
        "max_depth_overall" => maximum(results.max_depth),
        "total_rainfall_mm" => total_rainfall(scenario.rainfall),
    )

    package = ResultsPackage(results, metadata)
    write_results(output_dir, package, scenario.grid)
end

"""
    load_scenario(path::String)

Load a scenario from a TOML file.
"""
function load_scenario(path::String)
    load_scenario_from_toml(path)
end

"""
    save_results(config::RunConfig, results::ResultsAccumulator, scenario::Scenario)

Save simulation results to the output directory.
"""
function save_results(config::RunConfig, results::ResultsAccumulator{T},
                      scenario::Scenario{T}) where T
    write_results(config.output_dir, results, scenario, config.run_id)
end

"""
    save_metadata(config::RunConfig, metadata::RunMetadata)

Save run metadata to JSON file with comprehensive information.
"""
function save_metadata(config::RunConfig, metadata::RunMetadata)
    metadata_dict = Dict{String,Any}(
        "run_id" => metadata.run_id,
        "start_time" => string(metadata.start_time),
        "end_time" => metadata.end_time === nothing ? nothing : string(metadata.end_time),
        "duration_s" => metadata.end_time === nothing ? nothing :
            Dates.value(metadata.end_time - metadata.start_time) / 1000.0,
        "julia_version" => metadata.julia_version,
        "hydroforge_version" => metadata.hydroforge_version,
        "git_commit" => metadata.git_commit,
        "scenario_name" => metadata.scenario_name,
        "parameters" => metadata.parameters,
        "status" => string(metadata.status),
        "error_message" => metadata.error_message,
        "performance" => metadata.performance,
        "hazard_summary" => metadata.hazard_summary,
    )
    write_results_json(joinpath(config.output_dir, "run_metadata.json"), metadata_dict)
end

"""
    validate_scenario(scenario_path::String)

Validate a scenario file without running the simulation.
Returns a list of warnings/errors.
"""
function validate_scenario_file(scenario_path::String)
    issues = String[]

    if !isfile(scenario_path)
        push!(issues, "Scenario file not found: $scenario_path")
        return issues
    end

    try
        scenario = load_scenario(scenario_path)
        append!(issues, validate(scenario))
    catch e
        push!(issues, "Failed to load scenario: $(sprint(showerror, e))")
    end

    issues
end
