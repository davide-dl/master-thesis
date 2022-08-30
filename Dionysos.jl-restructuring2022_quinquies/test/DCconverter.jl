include("../src/Dionysos.jl")

import Plots

using StaticArrays
using ..Dionysos
using Profile
# using BenchmarkTools
DI = Dionysos
UT = DI.Utils
DO = DI.Domain
ST = DI.System
DM = DI.Map #domain_map
SR = DI.Relation #systemrelation
C = DI.Control

vs = 1.0; rL = 0.05; xL = 3.0; rC = 0.005; xC = 70.0; r0 = 1.0;

b = SVector(vs/xL, 0.0)
A1 = SMatrix{2,2}(-rL/xL, 0.0, 0.0, -1.0/xC/(r0+rC))
A2 = SMatrix{2,2}(-(rL+r0*rC/(r0+rC))/xL, 5.0*r0/(r0+rC)/xC,
    -r0/(r0+rC)/xL/5.0, -1.0/xC/(r0+rC))
F_sys = let b = b, A1 = A1, A2 = A2
    (x::SVector, u::SVector) -> u[1] == 1.0 ? A1*x + b : A2*x + b
end
tstep = 0.5
nsys = 5
sys1 = ST.NewSimpleSystemRK4(tstep,F_sys,nsys,SVector(0.0,0.0),SVector(0.0))

target = SVector(1.35,5.65)
accuracy = SVector(0.2,0.2)
safeSet = UT.HyperRectangle(target-accuracy, target+accuracy)
hx = SVector(5e-3, 5e-3)
# Xlims = UT.HyperRectangle(SVector(1.1, 5.4),SVector(1.5, 5.9))
Xlims = UT.HyperRectangle(SVector(1.0, 5.3),SVector(1.7, 6.0))
# Xlims = UT.HyperRectangle(SVector(1.15, 5.45), SVector(1.55, 5.85))
x0 = SVector(1.2, 5.6)
init = UT.HyperRectangle(SVector(1.195, 5.595), SVector(1.205, 5.605))

Xdom = DO.buildUniformDomain(Xlims,hx)
DO.addSet!(Xdom,safeSet)

# fig = Plots.plot(aspect_ratio = 1,legend = false)
# DO.plot_domain(Xdom,init=init,target=safeSet,fig=fig)

Xmap, Sdom = DM.buildNested2SymbolicMap(Xdom)
# println(length(DO.enum_elems(Sdom)))

Udom = DO.build1DListDomain([SVector{1,Float64}(0.0),SVector{1,Float64}(1.0)])
Umap, Adom = DM.buildStaticSymbolicMap(Udom)

sys2 = ST.buildEmptyDictProbAutomaton()
sys_rel = SR.buildConcrete2ProbSymb_Relation(Xmap,Umap,sys1,sys2,2)

# Profile.clear()
@time SR.buildAll!(sys_rel,simulate_borders=false,look_at_corners=false)
# Juno.profiler()
# println(sys_rel.sys2.nodes[319400])

@time best_actions, rewards = C.value_iteration_safety(sys_rel,target,accuracy)

C.simulateDC(sys_rel,x0,best_actions)
return
