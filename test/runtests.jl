using Test
using Aqua
using SIRcontrol

#
@testset verbose = true showtiming = true "SIRcontrol tests" begin
    for name in (:aqua, :default)
        @testset "$(name)" begin
            test_name = Symbol(:test_, name)
            include("$(test_name).jl")
            @eval $test_name()
        end
    end
end
