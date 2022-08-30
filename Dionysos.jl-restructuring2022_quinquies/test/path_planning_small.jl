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

function f(x, u)
    α = atan(tan(u[2])/2)
    return SVector{3}(
        u[1]*cos(α + x[3])/cos(α),
        u[1]*sin(α + x[3])/cos(α),
        u[1]*tan(u[2]))
end
tstep = 0.3
nsys = 5

sys1 = ST.NewSimpleSystemRK4(tstep,f,nsys,SVector(0.0,0.0,0.0),SVector(0.0,0.0))

angle_step = pi/5
Xlims = UT.HyperRectangle(SVector(8.2, 0.0,-pi), SVector(10.0,10.0,pi))
# Xlims = UT.HyperRectangle(SVector(8.5, 0.0,-pi), SVector(10.0,3.0,pi))
hx = SVector(0.4, 0.4, angle_step)
periods = SVector(SVector(0.0,0.0),SVector(0.0,0.0),SVector(-pi,2*pi))

Xdom = DO.buildPeriodicUniformDomain(Xlims,hx,periods)

X1_lb = [1.0, 2.2, 2.2, 3.4, 4.6, 5.8, 5.8, 7.0, 8.2, 8.4, 9.3, 8.4, 9.3, 8.4, 9.3];
X1_ub = [1.2, 2.4, 2.4, 3.6, 4.8, 6.0, 6.0, 7.2, 8.4, 9.3, 10.0, 9.3, 10.0, 9.3, 10.0];
X2_lb = [0.0, 0.0, 6.0, 0.0, 1.0, 0.0, 7.0, 1.0, 0.0, 8.2, 7.0, 5.8, 4.6, 3.4, 2.2];
X2_ub = [9.0, 5.0, 10.0, 9.0, 10.0, 6.0, 10.0, 10.0, 8.5, 8.6, 7.4, 6.2, 5.0, 3.8, 2.6];
for (x1lb, x2lb, x1ub, x2ub) in zip(X1_lb, X2_lb, X1_ub, X2_ub)
    box = UT.HyperRectangle(SVector(x1lb, x2lb, Xlims.lb[3]), SVector(x1ub, x2ub, Xlims.ub[3]))
    DO.removeObstacle!(Xdom,box)
end

I = UT.HyperRectangle(SVector(8.9, 8.9, -10), SVector(9.1, 9.1, 10))
# I = UT.HyperRectangle(SVector(8.5, 2.8, -pi-0.4), SVector(8.7, 3.0, pi+0.4))
# DO.addSet!(Xdom,I)
T = UT.HyperRectangle(SVector(9.2, 0.0,-10), SVector(10.0, 0.8,10))
DO.addSet!(Xdom,T)
#
Xmap, Sdom = DM.buildNested2SymbolicMap(Xdom)
#
fig = Plots.plot(aspect_ratio = 1,legend = false)
DO.plot_domain(Xdom,init=I,target=T,fig=fig)

U = UT.HyperRectangle(SVector(-1.0, -1.0), SVector(1.0, 1.0))

u0 = SVector(-1.0, 0.0)
h = SVector(2.0/5, 2.0/4)
Ugrid = DO.GridFree(u0, h)
Udom = DO.DomainList(Ugrid)
DO.add_set!(Udom, U, DO.OUTER)

Umap, Adom = DM.buildStaticSymbolicMap(Udom)

sys2 = ST.buildEmptyDictProbAutomaton()

sys_rel = SR.buildConcrete2ProbSymb_Relation(Xmap,Umap,sys1,sys2,1)

# Profile.clear()
# @time @profile SR.buildAll!(sys_rel,simulate_borders=false,look_at_corners=true)
# Juno.profiler()
@time SR.buildAll!(sys_rel,simulate_borders=false,look_at_corners=true)

T2 = DM.map1(sys_rel.xmap, T)

# @time best_actions, costs = C.backward_induction!(sys_rel,T2,lambda=1.0,MAXIT=500)
@time best_actions, costs = C.policy_iteration(sys_rel,T2,lambda=0.99)
# @time best_actions, costs = C.fast_bi_target(sys_rel,T)

# @time C.plot_costs!(sys_rel, costs,init=I,target=T,fig=fig)

xi = SVector(9.0, 9.0, 0.0)
C.simulate(sys_rel,xi,T,best_actions,costs,fig)

return
