include("../src/Dionysos.jl")

import Plots
import FileIO
using JLD2

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


tstep = 0.02
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

x = SVector(-1.0,-20.0,0.0,0.0)
u = SVector(-1.0)
for i=1:10
    global x = ST.computeTransition(sys1,x,u)
    # x_prev = x
    println(x)
end

return
