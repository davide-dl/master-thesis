module Control

using StaticArrays

using ..Utils
UT = Utils
using ..Domain
DO = Domain
using ..Map
DM = Map
using ..System
ST = System
using ..Relation
SR = Relation

@enum INCL_MODE INNER OUTER

include("MDP_control.jl")

end
