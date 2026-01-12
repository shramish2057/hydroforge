# Contributing

See [CONTRIBUTING.md](https://github.com/hydroforge/HydroForge/blob/main/CONTRIBUTING.md) for full guidelines.

## Quick Start

```bash
# Fork and clone
git clone https://github.com/shramish2057/HydroForge.git
cd HydroForge

# Install dependencies
julia --project -e 'using Pkg; Pkg.instantiate()'

# Run tests
julia --project -e 'using Pkg; Pkg.test()'

# Run demo
julia --project -e 'using HydroForge; HydroForge.run_demo()'
```

## Development Workflow

1. Create feature branch: `git checkout -b feature/my-feature`
2. Make changes
3. Run tests: `julia --project -e 'using Pkg; Pkg.test()'`
4. Format code: `julia -e 'using JuliaFormatter; format("src"); format("test")'`
5. Commit and push
6. Open pull request

## Code Style

We follow [Blue Style](https://github.com/invenia/BlueStyle):

- 4-space indentation
- 92-character line limit
- `snake_case` for functions
- `PascalCase` for types
- Docstrings for all public functions

## Testing

```julia
# Run all tests
using Pkg; Pkg.test()

# Run specific test file
include("test/test_types.jl")

# Run with coverage
Pkg.test(; coverage=true)
```

## Documentation

```julia
# Build docs locally
cd docs
julia --project make.jl

# Serve locally
julia -e 'using LiveServer; serve(dir="docs/build")'
```
