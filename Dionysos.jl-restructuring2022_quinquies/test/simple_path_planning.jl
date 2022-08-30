include("../src/Dionysos.jl")

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

tstep = 0.3

function F(x, u)
    return SVector{2}(
        x[1]+u[1],
        x[2]+u[2])
end

sys1 = ST.NewSimpleDiscreteSystem(tstep,F,SVector(0.0,0.0),SVector(0.0,0.0))

Xlims = UT.HyperRectangle(SVector(0.0, 0.0), SVector(10.0,10.0))
hx = SVector(0.5, 0.5)

Xdom = DO.buildUniformDomain(Xlims,hx)

X1_lb = [1.0, 2.2, 2.2, 3.4, 4.6, 5.8, 5.8, 7.0, 8.2, 8.4, 9.3, 8.4, 9.3, 8.4, 9.3];
X1_ub = [1.2, 2.4, 2.4, 3.6, 4.8, 6.0, 6.0, 7.2, 8.4, 9.3, 10.0, 9.3, 10.0, 9.3, 10.0];
X2_lb = [0.0, 0.0, 6.0, 0.0, 1.0, 0.0, 7.0, 1.0, 0.0, 8.2, 7.0, 5.8, 4.6, 3.4, 2.2];
X2_ub = [9.0, 5.0, 10.0, 9.0, 10.0, 6.0, 10.0, 10.0, 8.5, 8.6, 7.4, 6.2, 5.0, 3.8, 2.6];
for (x1lb, x2lb, x1ub, x2ub) in zip(X1_lb, X2_lb, X1_ub, X2_ub)
    box = UT.HyperRectangle(SVector(x1lb, x2lb), SVector(x1ub, x2ub))
    DO.removeObstacle!(Xdom,box)
end

I = UT.HyperRectangle(SVector(0.4, 0.4), SVector(0.6, 0.6))
DO.addSet!(Xdom, I)
T = UT.HyperRectangle(SVector(9.0, 0.0), SVector(10.0, 1.0))
DO.addSet!(Xdom,T)

Xmap, Sdom = DM.buildNested2SymbolicMap(Xdom)

fig = Plots.plot(aspect_ratio = 1,legend = false)
# DO.plot_domain(Xdom,init=I,target=T,fig=fig)

U = UT.HyperRectangle(SVector(-0.2, -0.2), SVector(0.2, 0.2))

u0 = SVector(0.0, 0.0)
h = SVector(0.19, 0.19)
Ugrid = DO.GridFree(u0, h)
Udom = DO.DomainList(Ugrid)
DO.add_set!(Udom, U, DO.OUTER)

Umap, Adom = DM.buildStaticSymbolicMap(Udom)

sys2 = ST.buildEmptyDictProbAutomaton()

sys_rel = SR.buildConcrete2ProbSymb_Relation(Xmap,Umap,sys1,sys2,8)

@time SR.buildFromInitSet!(sys_rel, I)

T2 = DM.map1(sys_rel.xmap, T)

xi = SVector(0.5, 0.5)
@time best_actions, values = C.value_iteration(sys_rel,T2,lambda=1.0,eps=0.001)

C.plot_values!(sys_rel, values,init=I,target=T,fig=fig)

function plot_actions!(r, actions)
    Plots.GR.setarrowsize(0.5)
    for c in DO.enum_centers(r.xmap.dom1[])
        s = DM.map1(r.xmap, c)
        if s!=0 && get(actions,s,0) != 0
            u = DM.map2(r.umap, actions[s])
            # line = Plots.Shape([(c[1], c[2]), (c[1]+u[1]/2, c[2]+u[2]/2)])
            # Plots.quiver!([c[1]],[c[2]],quiver=([u[1]],[u[2]]),color=:black,arrow=arrow(0.01, 0.01))
            # Plots.arrow!([c[1],c[2]],[u[1]*0.7,u[2]*0.7]))
            Plots.plot!([[c[1], c[1]+u[1]]], [[c[2], c[2]+u[2]]], arrow=true, opacity=0.5, color=:black)
        end
    end
    display(fig)
end
plot_actions!(sys_rel,best_actions)

C.simulate(sys_rel,xi,T,best_actions,fig)

return
