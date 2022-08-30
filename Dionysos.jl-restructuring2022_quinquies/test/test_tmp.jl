include("../src/Dionysos.jl")

using Profile
import Plots

using StaticArrays
using ..Dionysos
# using BenchmarkTools
DI = Dionysos
UT = DI.Utils
DO = DI.Domain
ST = DI.System
DM = DI.Map #domain_map
SR = DI.Relation #systemrelation
C = DI.Control

h = SVector(1,2,3)
x = SVector(1,4,10)

min(h./x...)
