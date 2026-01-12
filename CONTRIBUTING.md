# Contributing to HydroForge

Thank you for your interest in contributing to HydroForge! This document provides guidelines and instructions for contributing.

## Table of Contents

- [Code of Conduct](#code-of-conduct)
- [Getting Started](#getting-started)
- [Development Setup](#development-setup)
- [Code Style Guidelines](#code-style-guidelines)
- [Testing Requirements](#testing-requirements)
- [Pull Request Process](#pull-request-process)
- [Issue Reporting](#issue-reporting)

## Code of Conduct

This project adheres to the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md). By participating, you are expected to uphold this code. Please report unacceptable behavior to the project maintainers.

## Getting Started

1. Fork the repository on GitHub
2. Clone your fork locally
3. Set up the development environment (see below)
4. Create a branch for your changes
5. Make your changes and commit them
6. Push to your fork and submit a pull request

## Development Setup

### Prerequisites

- Julia 1.10 or higher
- Git

### Installation

```bash
# Clone your fork
git clone https://github.com/shramish2057/HydroForge.git
cd HydroForge

# Install dependencies
julia --project -e 'using Pkg; Pkg.instantiate()'

# Run tests to verify setup
julia --project -e 'using Pkg; Pkg.test()'
```

### Running the Demo

```julia
using HydroForge
HydroForge.run_demo()
```

## Code Style Guidelines

We follow the [Blue Style](https://github.com/invenia/BlueStyle) for Julia code formatting.

### Formatting

Use JuliaFormatter to format your code:

```julia
using JuliaFormatter
format("src")
format("test")
```

### Key Style Points

- **Indentation**: 4 spaces (no tabs)
- **Line length**: Maximum 92 characters
- **Naming conventions**:
  - `snake_case` for functions and variables
  - `PascalCase` for types and modules
  - `SCREAMING_SNAKE_CASE` for constants
- **Documentation**: All public functions must have docstrings
- **Type annotations**: Use for function signatures, optional for local variables

### Example

```julia
"""
    compute_flux(depth::T, velocity::T, grid::Grid{T}) where T<:AbstractFloat

Compute the discharge flux at cell faces.

# Arguments
- `depth`: Water depth (m)
- `velocity`: Flow velocity (m/s)
- `grid`: Computational grid

# Returns
- Computed flux value (m²/s)
"""
function compute_flux(depth::T, velocity::T, grid::Grid{T}) where T<:AbstractFloat
    return depth * velocity * grid.dx
end
```

## Testing Requirements

All contributions must include appropriate tests.

### Running Tests

```bash
# Run all tests
julia --project -e 'using Pkg; Pkg.test()'

# Run specific test file
julia --project test/test_types.jl
```

### Test Guidelines

1. **Unit tests**: Test individual functions in isolation
2. **Integration tests**: Test component interactions
3. **Edge cases**: Include tests for boundary conditions
4. **Coverage**: Aim for >80% code coverage for new code

### Test Structure

```julia
@testset "Grid Tests" begin
    @testset "Constructor" begin
        grid = Grid(100, 100, 10.0, 10.0)
        @test grid.nx == 100
        @test grid.ny == 100
    end

    @testset "Edge Cases" begin
        @test_throws ArgumentError Grid(-1, 100, 10.0, 10.0)
    end
end
```

## Pull Request Process

### Before Submitting

1. Ensure all tests pass locally
2. Format your code with JuliaFormatter
3. Update documentation if needed
4. Add yourself to CONTRIBUTORS.md (if applicable)

### PR Guidelines

1. **Title**: Use a clear, descriptive title
2. **Description**: Explain what changes you made and why
3. **Issue reference**: Link to related issues using `Fixes #123` or `Relates to #456`
4. **Small PRs**: Keep pull requests focused and reasonably sized
5. **Draft PRs**: Use draft PRs for work-in-progress

### Review Process

1. A maintainer will review your PR
2. Address any requested changes
3. Once approved, a maintainer will merge your PR

## Issue Reporting

### Bug Reports

When reporting bugs, please include:

1. **Description**: Clear description of the bug
2. **Reproduction steps**: Minimal code to reproduce the issue
3. **Expected behavior**: What you expected to happen
4. **Actual behavior**: What actually happened
5. **Environment**: Julia version, OS, HydroForge version
6. **Error messages**: Full error output if applicable

### Feature Requests

When requesting features, please include:

1. **Use case**: Why this feature would be useful
2. **Proposed solution**: How you envision it working
3. **Alternatives**: Any alternative solutions you've considered

### Questions

For questions about the code or usage:

1. Check existing documentation and issues first
2. Use GitHub Discussions for general questions
3. Open an issue for specific technical questions

## Getting Help

- **Documentation**: See the `docs/` directory
- **Issues**: Search existing issues or create a new one
- **Discussions**: Use GitHub Discussions for general questions

Thank you for contributing to HydroForge!
