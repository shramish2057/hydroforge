# HydroForge Mass Balance Module
# Mass balance tracking and error computation

"""
    MassBalance{T<:AbstractFloat}

Tracks mass balance components throughout simulation.

# Fields
- `initial_volume::T`: Initial water volume in domain (m³)
- `rainfall_volume::T`: Cumulative rainfall added (m³)
- `outflow_volume::T`: Cumulative outflow through boundaries (m³)
- `infiltration_volume::T`: Cumulative infiltration losses (m³)
- `current_volume::T`: Current water volume in domain (m³)
"""
mutable struct MassBalance{T<:AbstractFloat}
    initial_volume::T
    rainfall_volume::T
    outflow_volume::T
    infiltration_volume::T
    current_volume::T
end

"""
    MassBalance(T=Float64)

Create a new mass balance tracker with zero initial values.
"""
function MassBalance(::Type{T}=Float64) where T<:AbstractFloat
    MassBalance{T}(zero(T), zero(T), zero(T), zero(T), zero(T))
end

"""
    MassBalance(state::SimulationState, grid::Grid)

Create a mass balance tracker initialized with current state volume.
"""
function MassBalance(state::SimulationState{T}, grid::Grid{T}) where T
    vol = total_volume(state, grid)
    MassBalance{T}(vol, zero(T), zero(T), zero(T), vol)
end

"""
    reset!(mb::MassBalance, state::SimulationState, grid::Grid)

Reset mass balance tracker with current state as initial condition.
"""
function reset!(mb::MassBalance{T}, state::SimulationState{T}, grid::Grid{T}) where T
    vol = total_volume(state, grid)
    mb.initial_volume = vol
    mb.rainfall_volume = zero(T)
    mb.outflow_volume = zero(T)
    mb.infiltration_volume = zero(T)
    mb.current_volume = vol
    mb
end

"""
    update_volume!(mb::MassBalance, state::SimulationState, grid::Grid)

Update current volume from simulation state.
"""
function update_volume!(mb::MassBalance{T}, state::SimulationState{T}, grid::Grid{T}) where T
    mb.current_volume = total_volume(state, grid)
    mb
end

"""
    add_rainfall!(mb::MassBalance, volume::Real)

Record rainfall volume added to domain.
"""
function add_rainfall!(mb::MassBalance{T}, volume::Real) where T
    mb.rainfall_volume += T(volume)
    mb
end

"""
    add_outflow!(mb::MassBalance, volume::Real)

Record outflow volume leaving domain.
"""
function add_outflow!(mb::MassBalance{T}, volume::Real) where T
    mb.outflow_volume += T(volume)
    mb
end

"""
    add_infiltration!(mb::MassBalance, volume::Real)

Record infiltration volume.
"""
function add_infiltration!(mb::MassBalance{T}, volume::Real) where T
    mb.infiltration_volume += T(volume)
    mb
end

"""
    mass_error(mb::MassBalance)

Compute mass balance error (m³).

Error = initial + rainfall - outflow - infiltration - current

Should be zero for perfect mass conservation.
"""
function mass_error(mb::MassBalance{T}) where T
    mb.initial_volume + mb.rainfall_volume - mb.outflow_volume -
        mb.infiltration_volume - mb.current_volume
end

"""
    relative_mass_error(mb::MassBalance)

Compute relative mass balance error (dimensionless).

Returns error / reference_volume, where reference is max of initial or total input.
Returns 0 if no water has entered the system.
"""
function relative_mass_error(mb::MassBalance{T}) where T
    error = mass_error(mb)

    # Reference volume: max of initial volume or total inputs
    reference = max(mb.initial_volume, mb.rainfall_volume)

    if reference < eps(T)
        return zero(T)
    end

    abs(error) / reference
end

"""
    compute_mass_balance(state::SimulationState, grid::Grid,
                         initial_volume::Real, rainfall_volume::Real,
                         outflow_volume::Real)

Compute mass balance error for a simulation.

# Returns
- `error`: Absolute mass error (m³)
- `relative_error`: Relative error (dimensionless)
"""
function compute_mass_balance(state::SimulationState{T}, grid::Grid{T},
                              initial_volume::Real, rainfall_volume::Real,
                              outflow_volume::Real) where T
    current = total_volume(state, grid)

    # Mass balance: initial + inputs = current + outputs
    # Error = initial + rainfall - outflow - current
    error = T(initial_volume) + T(rainfall_volume) - T(outflow_volume) - current

    # Relative error
    reference = max(T(initial_volume), T(rainfall_volume))
    if reference < eps(T)
        relative = zero(T)
    else
        relative = abs(error) / reference
    end

    (error=error, relative_error=relative)
end

"""
    check_mass_balance(mb::MassBalance; tolerance=0.01)

Check if mass balance error is within tolerance.

# Arguments
- `mb`: Mass balance tracker
- `tolerance`: Maximum acceptable relative error (default 1%)

# Returns
- `Bool`: true if error is within tolerance
"""
function check_mass_balance(mb::MassBalance; tolerance::Real=0.01)
    relative_mass_error(mb) <= tolerance
end
