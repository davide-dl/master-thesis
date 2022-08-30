include("../src/Dionysos.jl")

LOAD_REL = false
LOAD_CTRL = false
SAVE_REL = false
SAVE_CTRL = false

# SAVE_REL = true
# LOAD_REL = true

# LOAD_CTRL = true
SAVE_CTRL = true

import Plots
import FileIO
using JLD2
using Profile
# using JLD

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


# tstep = 0.03 #CHANGED
tstep = 0.04
function F(state, u)
    gravity = 9.8
    masscart = 1.0
    masspole = 0.1
    total_mass = masspole + masscart
    pole_length = 0.5  # actually half the pole's length
    polemass_length = masspole * pole_length
    force_mag = 10.0
    tau = tstep
    force = force_mag*u[1]
    x = state[1]
    theta = state[2]
    x_dot = state[3]
    theta_dot = state[4]

    temp = (force + polemass_length * (theta_dot *pi/180)^2 * sin(theta *pi/180)) / total_mass
    thetaacc = (gravity * sin(theta *pi/180) - cos(theta *pi/180) * temp) / (pole_length * (4.0 / 3.0 - masspole * cos(theta *pi/180)^2 / total_mass)) * 180/pi
    xacc = temp - polemass_length * thetaacc *pi/180 * cos(theta *pi/180) / total_mass
    return SVector(
        x + tau * x_dot,
        theta + tau * theta_dot,
        x_dot + tau * xacc,
        theta_dot + tau * thetaacc)
end

sys1 = ST.NewSimpleDiscreteSystem(tstep,F,SVector(0.0,0.0,0.0,0.0),SVector(0.0))

x_th = 5.2
theta_th = 45.0
x_dot_th = 8.0
theta_dot_th = 4.0 * 180/pi
th_left = SVector(x_th,12.0,x_dot_th,theta_dot_th)
th_right = SVector(0.8,theta_th,x_dot_th,theta_dot_th)
th = SVector(x_th,theta_th,x_dot_th,theta_dot_th)
th_easy = SVector(0.9,18.0,x_dot_th,theta_dot_th)

x_step = 0.4
theta_step = 3.0
x_dot_step = x_dot_th/15
theta_dot_step = theta_dot_th/15
step = SVector(x_step,theta_step,x_dot_step,theta_dot_step)
step_easy = SVector(0.3,3.0,x_dot_step,theta_dot_step)

Xlims = UT.HyperRectangle(-th_left, th_right) #NOTE
Xdom = DO.buildUniformDomain(Xlims,step) #NOTE

# I = UT.HyperRectangle(SVector(8.5, 2.8, -pi-0.4), SVector(8.7, 3.0, pi+0.4))
# DO.addSet!(Xdom,I)
target = SVector(x_step,theta_step,x_dot_th,theta_dot_th)
target_easy = SVector(0.3,3.0,x_dot_th,theta_dot_th)
T = UT.HyperRectangle(-target,target) #NOTE
DO.addSet!(Xdom,T)

fig = Plots.plot(aspect_ratio = step_easy[1]/step_easy[2],legend = false)
# fig = Plots.plot(aspect_ratio=x_step/theta_step,
#                 xlims = (-4,1),
#                 xticks = -4:1.0:1,
#                 yticks = -10:10.0:30,
#                 ylims = (-12,30),
#                 legend = false)
DO.plot_domain(Xdom,target=T,fig=fig)

Xmap, Sdom = DM.buildNested2SymbolicMap(Xdom)

Udom = DO.build1DListDomain([SVector(-1.0),SVector(1.0)])

Umap, Adom = DM.buildStaticSymbolicMap(Udom)

sys2 = ST.buildEmptyDictProbAutomaton()

sys_rel = SR.buildConcrete2ProbSymb_Relation(Xmap,Umap,sys1,sys2,1)

if LOAD_REL
    loaded_rel = FileIO.load("relation.jld2")
    sys_rel = loaded_rel["r"]
else
    @time SR.buildAll!(sys_rel,simulate_borders=false,look_at_corners=false,var_tstep=0)
    # @time SR.buildAll_easy!(sys_rel)
    # if SAVE_REL
    #     FileIO.save("relation.jld2", Dict("r"=>sys_rel))
    #     # save_object("relation.jld2", sys_rel)
    #     # @save "relation.jld2" sys_rel
    #     println("relation saved")
    # end
end

# s = DM.map1(sys_rel.xmap,SVector(0.0,0.0,pi/24,0.0))
# node = sys_rel.sys2.nodes[s]
# action = node.actions[1]
# println(action)

x_bound = 2.2
theta_bound = 12.0
bounds = SVector(x_bound,theta_bound,x_dot_th,theta_dot_th)
bounds_easy = SVector(0.6,12.0,x_dot_th,theta_dot_th)
B = UT.HyperRectangle(-th_left,th_right) #NOTE
T2 = DM.map1(sys_rel.xmap, T)

# iterations=80 #CHANGED
iterations=100
discount=1.0
if LOAD_CTRL
    loaded_ctrl = FileIO.load("controller.jld2")
    best_actions = loaded_ctrl["c"]
    rewards = loaded_ctrl["v"]
else
    # Profile.clear()
    # @time @profile best_actions, rewards = C.vi_target_cart(sys_rel,T,B,MAXIT=iterations)
    # Juno.profiler()
    # @time best_actions, rewards = C.vi_target_cart(sys_rel,T,B,lambda=discount,MAXIT=iterations)
    # @time best_actions, rewards = C.bi_target(sys_rel,T,MAXIT=iterations)
    @time best_actions, costs = C.backward_induction!(sys_rel,T2,lambda=discount,MAXIT=iterations,var_tstep=false)
    # @time best_actions, costs = C.fast_bi_target(sys_rel,T)
    rewards = costs
    if SAVE_CTRL
        FileIO.save("controller.jld2", Dict("c"=>best_actions, "v"=>rewards))
        println("controller saved")
    end
end

# @time C.plot_values!(sys_rel,rewards,target=T,fig=fig,scale=iterations*discount)
# @time C.plot_costs!(sys_rel,rewards,target=T,fig=fig)

xi = SVector(-3.0,-10.0,0.0,0.0)
# si = DM.map1(sys_rel.xmap,xi)
# actions = sys_rel.sys2.nodes[si].actions
# for (a, action) in actions
#     println("action ", a)
#     for (s,t) in action.trans
#         println("state ",s, " prob ", t.prob, " x ", DM.map2(sys_rel.xmap,s))
#     end
# end

# C.simulate_cartpole(sys_rel,xi,T,best_actions,rewards,fig)
C.simulate_cartpole(sys_rel,xi,T,best_actions,costs,fig)

return
