module Map

using StaticArrays
using ..Domain
DO = Domain
using ..Utils
UT = Utils

abstract type AbstractDomainMap{N1,T1,N2,T2} end

@enum INCL_MODE INNER OUTER

include("domain_map.jl")

end
