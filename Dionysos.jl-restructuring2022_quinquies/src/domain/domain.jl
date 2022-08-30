module Domain

# using LinearAlgebra
# using ProgressMeter
using StaticArrays
# using Base.Cartesian
# using HybridSystems
using ..Utils
UT = Utils

@enum INCL_MODE INNER OUTER

"""
    DomainType{N,T}

Description
"""
abstract type DomainType{N,T} end

include("grid.jl")
include("domain_list.jl")
include("NestedDomain/nested_domain.jl") #Davide
include("list_domain.jl") #Davide
include("SimpleSymbolicDomain/simple_symbolic_domain.jl") #Davide
end  # module Domain
