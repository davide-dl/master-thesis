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

function F(state, u, tstep)
    gravity = 9.8
    pole_length = 0.5  # actually half the pole's length
    acc = 20.0

    theta = state[1]
    theta_dot = state[2]

    applied_acc = (pole_length * acc * u[1] * cos(theta*pi/180)) * 180/pi
    g_acc = (pole_length * gravity * sin(theta*pi/180)) * 180/pi
    theta_acc = applied_acc + g_acc

    return SVector(
        theta + tstep * theta_dot,
        theta_dot + tstep * theta_acc)
end

function f(state, u)
    gravity = 9.8
    pole_length = 0.5  # actually half the pole's length
    acc = 20.0

    theta = state[1]
    theta_dot = state[2]

    applied_acc = (pole_length * acc * u[1] * cos(theta*pi/180)) * 180/pi
    g_acc = (pole_length * gravity * sin(theta*pi/180)) * 180/pi
    theta_acc = applied_acc + g_acc

    return SVector(
        theta_dot,
        theta_acc)
end

nsys = 5
min_tstep = 0.04
max_tstep = 0.04

RKsys = ST.NewVarTSystemRK4(f,nsys,SVector(0.0,0.0),SVector(0.0),min_tstep,max_tstep)
# sys1 = ST.NewVarTDiscreteSystem(f,F,SVector(0.0,0.0),SVector(0.0),min_tstep,max_tstep)
sys1 = RKsys

theta_th = 45.0
theta_dot_th = 4.0 * 180/pi
th = SVector(theta_th,theta_dot_th)

theta_step = 3.0
theta_dot_step = 2*theta_dot_th/31
step = SVector(theta_step,theta_dot_step)

Xlims = UT.HyperRectangle(-th, th)
Xdom = DO.buildUniformDomain(Xlims,step) #NOTE

target = SVector(theta_step,theta_dot_th)
T = UT.HyperRectangle(-target,target)
DO.addSet!(Xdom,T)

fig = Plots.plot(aspect_ratio = theta_step/theta_dot_step, legend = false)
# DO.plot_domain(Xdom,target=T,fig=fig)

Xmap, Sdom = DM.buildNested2SymbolicMap(Xdom)

Udom = DO.build1DListDomain([SVector(-1.0),SVector(1.0)])
Umap, Adom = DM.buildStaticSymbolicMap(Udom)

sys2 = ST.buildEmptyDictProbAutomaton()

sys_rel = SR.buildConcrete2ProbSymb_Relation(Xmap,Umap,sys1,sys2,10)

@time SR.buildAll!(sys_rel,simulate_borders=true,look_at_corners=false,var_tstep=1)

T2 = DM.map1(sys_rel.xmap, T)

iterations=10
discount=1.0
@time best_actions, costs = C.backward_induction!(sys_rel,T2,lambda=discount,MAXIT=iterations,var_tstep=true)

@time C.plot_costs!(sys_rel,costs,fig=fig)

xi = SVector(-20.0,0.0)

C.simulateRK4(RKsys, sys_rel,xi,T,best_actions,costs,fig,tstep=0,render_each_step=true)

xw = SVector(8.0,-210.0)
sw = DM.map1(sys_rel.xmap,xw)
println(sw)
println("cost ", costs[sw])
for (st,t) in sys_rel.sys2.nodes[sw].actions[2].trans
    xt = DM.map2(sys_rel.xmap,st)
    println("new_x = ", xt)
    println("new_s = ", DM.map1(sys_rel.xmap,xt))
    println("cost ", costs[st])
end

return
