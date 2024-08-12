using VertexPostprocessing
using Test

@testset "ChannelTransforms" begin
    include("ChannelTransforms.jl")
end

@testset "Expand test" begin
    include("expand_test.jl")
end

@testset "GFTools" begin
    include("GFTools.jl")
end

@testset "FullRun" begin
    empty!(ARGS)
    push!(ARGS, abspath("./test_data/data_s1/grid_b5_f5_s1.jld2"))
    push!(ARGS, abspath("./test_data/data_s1"))
    include("../scripts/expand_vertex.jl")
end
