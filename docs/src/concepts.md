# Concepts

## Physical Model

HydroForge simulates **pluvial (rainfall-driven) urban flooding** using the 2D shallow water equations with the local inertial approximation.

### Governing Equations

The shallow water equations describe conservation of mass and momentum for free surface flows:

**Mass Conservation (Continuity):**
```math
\frac{\partial h}{\partial t} + \frac{\partial q_x}{\partial x} + \frac{\partial q_y}{\partial y} = R - I
```

**Momentum (Local Inertial Form):**
```math
\frac{\partial q_x}{\partial t} = -gh\frac{\partial \eta}{\partial x} - \frac{gn^2|q_x|q_x}{h^{10/3}}
```

```math
\frac{\partial q_y}{\partial t} = -gh\frac{\partial \eta}{\partial y} - \frac{gn^2|q_y|q_y}{h^{10/3}}
```

### Key Variables

| Variable | Symbol | Units | Description |
|----------|--------|-------|-------------|
| Water depth | ``h`` | m | Depth of water above ground |
| Unit discharge (x) | ``q_x`` | m²/s | Flow per unit width in x |
| Unit discharge (y) | ``q_y`` | m²/s | Flow per unit width in y |
| Water surface elevation | ``\eta = h + z`` | m | Water level above datum |
| Bed elevation | ``z`` | m | Ground elevation |
| Rainfall rate | ``R`` | m/s | Rainfall intensity |
| Infiltration rate | ``I`` | m/s | Water loss to ground |
| Manning's n | ``n`` | s/m^(1/3) | Surface roughness |

### Why Local Inertial?

The **local inertial approximation** neglects the advection terms in the full shallow water equations. This provides:

1. **Numerical stability**: No need for Riemann solvers
2. **Efficiency**: Simple explicit time-stepping
3. **Accuracy**: Suitable for subcritical urban flows
4. **Robustness**: Handles wetting/drying well

## Numerical Method

### Spatial Discretization

HydroForge uses a **structured regular grid** with:

- Cell-centered depth and discharge storage
- Face-based flux computation
- Central differences for gradients

### Temporal Discretization

The solver uses **semi-implicit time-stepping**:

1. Fluxes are computed explicitly
2. Friction is treated semi-implicitly for stability
3. Adaptive timestep based on CFL condition

### CFL Condition

The timestep is limited by:

```math
\Delta t \leq C_{CFL} \cdot \frac{\Delta x}{\sqrt{gh} + |u|}
```

where ``C_{CFL} \approx 0.7`` is the CFL number.

### Wetting and Drying

Cells transition between wet and dry states:

- **Wet cell**: ``h > h_{min}`` (default 1 mm)
- **Dry cell**: ``h \leq h_{min}``
- Fluxes are zero at dry cell faces

## Boundary Conditions

| Type | Description | Use Case |
|------|-------------|----------|
| Closed | No flow across boundary | Domain edges, walls |
| Open | Free outflow | Natural drainage outlets |
| Fixed depth | Prescribed water level | Tidal boundaries |

## Manning's Roughness

Typical values for urban surfaces:

| Surface Type | Manning's n |
|--------------|-------------|
| Smooth concrete | 0.011-0.013 |
| Asphalt | 0.013-0.016 |
| Gravel | 0.020-0.030 |
| Short grass | 0.030-0.035 |
| Dense vegetation | 0.040-0.080 |
| Buildings (blockage) | 0.100-0.200 |

## Limitations

Current limitations of HydroForge:

1. **2D only**: No 1D drainage network coupling (planned)
2. **Single-phase**: Water only, no sediment transport
3. **Uniform rainfall**: Spatially uniform precipitation (spatial rainfall planned)
4. **No infiltration**: Green-Ampt model planned for future
5. **Subcritical flows**: Not validated for supercritical conditions
