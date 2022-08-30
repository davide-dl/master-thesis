module System

using StaticArrays
using HybridSystems

abstract type AbstractSystem{Nx,Nu,T} end

include("controlsystem.jl")
include("SimpleSystem/simple_system.jl")
include("ProbabilisticSymbolicSystem/prob_symbolic_system.jl")

end
