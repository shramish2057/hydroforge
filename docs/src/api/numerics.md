# Numerics API

Numerical methods for shallow water computation.

## Timestep Control

```@docs
compute_dt
compute_dt_array
check_cfl
```

## Flux Computation

```@docs
compute_velocity
water_surface_elevation
water_surface_elevation!
face_depth_x
face_depth_y
compute_flux_x!
compute_flux_y!
```

## Boundary Conditions

```@docs
BoundaryType
apply_boundaries!
apply_closed_boundaries!
apply_open_boundaries!
enforce_positive_depth!
```

## Solver

```@docs
SimulationWorkspace
step!
update_depth!
run_simulation!
```
