include("../../src/Dionysos.jl")

using StaticArrays
using ..Dionysos
DI = Dionysos
UT = DI.Utils
DO = DI.Domain
# const ST = DI.System
# const CO = DI.Control
# const SY = DI.Symbolic

a = SVector(0.0,0.0)
b = SVector(10.0,10.0)
lims = UT.HyperRectangle(a,b)
h = SVector(1.0,1.0)

dom = DO.buildUniformDomain(lims,h)
DO.removeObstacle!(dom,UT.HyperRectangle(SVector(0.5,0.5),SVector(3.0,3.0)))
DO.removeObstacle!(dom,UT.HyperRectangle(SVector(0.7,0.5),SVector(7.5,3.0)))

DO.refineCell!(dom,SVector(5,5),SVector(0.1,0.229))

A = UT.HyperRectangle(SVector(5.1,5.2),SVector(7.0,8.0))
DO.addSet!(dom,A)
B = UT.HyperRectangle(SVector(3.5,7.5),SVector(9.5,10.0))
DO.addSet!(dom,B)

# for (key, value) in dom.subdomains
#     println(key)
# end

DO.plot_domain(dom)
