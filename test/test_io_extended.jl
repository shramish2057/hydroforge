# Extended IO Tests for HydroForge
# Tests for GeoJSON reader, Network TOML loader, and ESRI ASCII Grid reader

using Test
using HydroForge

@testset "Extended IO Tests" begin
    fixtures_dir = joinpath(@__DIR__, "fixtures")

    @testset "GeoJSON Reader" begin
        geojson_path = joinpath(fixtures_dir, "test_points.geojson")

        @testset "Read points as coordinates" begin
            points = read_points_geojson(geojson_path)

            # Should have 4 points total (2 Point + 2 from MultiPoint)
            @test length(points) == 4

            # Check first point
            @test points[1] == (100.5, 50.5)

            # Check second point
            @test points[2] == (200.5, 150.5)

            # Check MultiPoint coordinates
            @test points[3] == (300.0, 200.0)
            @test points[4] == (400.0, 250.0)
        end

        @testset "Read points with grid conversion" begin
            grid = Grid(50, 50, 10.0)  # 50x50 grid with 10m cells

            cell_indices = read_points_geojson(geojson_path, grid)

            # Should have 4 cell indices
            @test length(cell_indices) == 4

            # Check conversion (100.5 / 10) + 1 = 11, (50.5 / 10) + 1 = 6
            @test cell_indices[1] == (11, 6)
        end

        @testset "point_to_cell conversion" begin
            grid = Grid(100, 100, 10.0, 10.0, 0.0, 0.0)

            # Point at origin - floor(5/10) + 1 = 1
            @test point_to_cell(5.0, 5.0, grid) == (1, 1)

            # Point at center - floor(500/10) + 1 = 51
            @test point_to_cell(500.0, 500.0, grid) == (51, 51)

            # Point clamped to bounds
            @test point_to_cell(-100.0, -100.0, grid) == (1, 1)
            @test point_to_cell(2000.0, 2000.0, grid) == (100, 100)
        end

        @testset "Error handling" begin
            @test_throws ArgumentError read_points_geojson("nonexistent.geojson")
        end
    end

    @testset "Network TOML Loader" begin
        network_path = joinpath(fixtures_dir, "test_network.toml")

        @testset "Load network from TOML" begin
            network = load_network_from_toml(network_path)

            # Check junctions
            @test n_junctions(network) == 3
            @test network.junctions[1].id == 1
            @test network.junctions[1].x == 100.0
            @test network.junctions[1].invert == -2.0
            @test network.junctions[1].junction_type == MANHOLE
            @test network.junctions[3].junction_type == OUTFALL

            # Check pipes
            @test n_pipes(network) == 2
            @test network.pipes[1].id == 1
            @test network.pipes[1].upstream_node == 1
            @test network.pipes[1].downstream_node == 2
            @test network.pipes[1].section isa CircularPipe
            @test network.pipes[1].section.diameter == 0.45

            @test network.pipes[2].section isa RectangularPipe
            @test network.pipes[2].section.width == 0.6
            @test network.pipes[2].section.height == 0.4

            # Check inlets
            @test n_inlets(network) == 2
            @test network.inlets[1].inlet_type == GRATE
            @test network.inlets[2].inlet_type == CURB
            @test network.inlets[2].length == 1.2

            # Check outlets
            @test length(network.outlets) == 1
            @test network.outlets[1].outlet_type == :FREE
        end

        @testset "Network validation" begin
            network = load_network_from_toml(network_path)
            issues = validate(network)

            # Network should be valid
            @test isempty(issues)
        end

        @testset "Read raw TOML" begin
            config = read_network_toml(network_path)

            @test haskey(config, "network")
            @test config["network"]["name"] == "Test Network"
            @test haskey(config, "junctions")
            @test haskey(config, "pipes")
        end

        @testset "Error handling" begin
            @test_throws ArgumentError load_network_from_toml("nonexistent.toml")
        end
    end

    @testset "ESRI ASCII Grid Reader" begin
        asc_path = joinpath(fixtures_dir, "test_dem.asc")

        @testset "Read ASCII grid" begin
            data, metadata = read_esri_ascii(asc_path)

            # Check dimensions
            @test size(data) == (10, 10)

            # Check metadata
            @test metadata["ncols"] == 10
            @test metadata["nrows"] == 10
            @test metadata["dx"] ≈ 10.0
            @test metadata["dy"] ≈ 10.0
            @test metadata["x_origin"] == 0.0
            @test metadata["y_origin"] == 0.0

            # Check data values
            @test data[1, 1] ≈ 0.0
            @test data[1, 10] ≈ 0.9
            @test data[10, 10] ≈ 1.8
        end

        @testset "read_geotiff with .asc extension" begin
            data, metadata = read_geotiff(asc_path)

            # Should work via read_geotiff
            @test size(data) == (10, 10)
            @test haskey(metadata, "dx")
        end

        @testset "Error handling" begin
            @test_throws ArgumentError read_esri_ascii("nonexistent.asc")
        end
    end

    @testset "read_geotiff format detection" begin
        # Create temp files to test format detection
        temp_dir = mktempdir()

        @testset "Unsupported format" begin
            # Create a dummy file with unsupported extension
            xyz_path = joinpath(temp_dir, "test.xyz")
            write(xyz_path, "dummy content")
            @test_throws ErrorException read_geotiff(xyz_path)
        end

        @testset "GeoTIFF without ArchGDAL" begin
            # Create a dummy .tif file to test ArchGDAL requirement
            tif_path = joinpath(temp_dir, "test.tif")
            write(tif_path, "dummy tiff content")
            # Should error about ArchGDAL not being available
            @test_throws ErrorException read_geotiff(tif_path)
        end

        @testset "File not found" begin
            @test_throws ArgumentError read_geotiff("nonexistent.asc")
        end

        # Cleanup
        rm(temp_dir, recursive=true)
    end
end
