module Relation

using StaticArrays

using ..Utils
UT = Utils
using ..Domain
DO = Domain
using ..Map
using ..System

abstract type AbstractRelation{Nx1,Nu1,T1,Nx2,Nu2,T2} end

@enum INCL_MODE INNER OUTER

include("concrete2probsymb.jl")

end
