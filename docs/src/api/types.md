# Types API

Core data types for HydroForge simulations.

## Grid

```@docs
Grid
cell_area
total_area
extent
cell_centers_x
cell_centers_y
cell_index
```

## Topography

```@docs
Topography
min_elevation
max_elevation
elevation_range
mean_roughness
```

## Simulation State

```@docs
SimulationState
total_volume
max_depth
min_depth
mean_depth
wet_cells
wet_fraction
max_velocity
copy_state
reset!
```

## Parameters

```@docs
SimulationParameters
validate
```

## Rainfall

```@docs
RainfallEvent
rainfall_rate
rainfall_rate_ms
total_rainfall
duration
peak_intensity
```

## Scenario

```@docs
Scenario
```

## Results

```@docs
ResultsPackage
ResultsAccumulator
```
