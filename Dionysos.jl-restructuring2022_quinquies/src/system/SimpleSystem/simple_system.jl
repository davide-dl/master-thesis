import Random
import Distributions

struct SimpleSystem{Nx,Nu,T} <: AbstractSystem{Nx,Nu,T}
    tstep::Float64
    F::Function #function x(t+1)=F(x(t),u(t))
    # F_inv::Function
end
struct VarTSystem{Nx,Nu,T} <: AbstractSystem{Nx,Nu,T}
    min_tstep::Float64
    max_tstep::Float64
    f::Function #function x'(t)=f(x(t),u(t))
    F::Function #function x(t+1)=F(x(t),u(t))
    # F_inv::Function
end

struct VarTDiscreteSystem{Nx,Nu,T} <: AbstractSystem{Nx,Nu,T}
    min_tstep::Float64
    max_tstep::Float64
    f::Function #function x'(t)=f(x(t),u(t))
    F::Function #function x(t+1)=F(x(t),u(t))
    # F_inv::Function
end

function NewVarTSystemRK4(f_sys::Function, nsys::Int, xtmp::SVector{Nx,T}, utmp::SVector{Nu,T}, min_tstep=0,max_tstep=Inf) where {Nx,Nu,T}
    sys_map = let nsys = nsys
        (x::SVector{Nx,T}, u::SVector{Nu,T}, tstep::Real) ->
        RungeKutta4(f_sys, x, u, tstep, nsys)::SVector{Nx,T} #use Vector instead?
    end
    # sys_inv_map = let nsys = nsys
    #     (x::SVector{N,T}, u, tstep) ->
    #         RungeKutta4(F_sys, x, u, -tstep, nsys)::SVector{N,T}
    # end
    return VarTSystem{Nx,Nu,T}(min_tstep,max_tstep,f_sys,sys_map)
end

function NewNoisyVarTSystemRK4(f_sys::Function, nsys::Int, xtmp::SVector{Nx,T}, utmp::SVector{Nu,T}, sigma, min_tstep=0,max_tstep=Inf) where {Nx,Nu,T}
    Random.seed!(1)
    noise = Vector{T}(undef,Nx)
    for i = 1 : Nx
        d = Distributions.Normal(0,sigma[i])
        noise[i] = rand(d)
    end
    sys_map = let nsys = nsys
        (x::SVector{Nx,T}, u::SVector{Nu,T}, tstep::Real) ->
        RungeKutta4(f_sys, x, u, tstep, nsys)::SVector{Nx,T}+SVector{Nx}(noise) #use Vector instead?
    end
    # sys_inv_map = let nsys = nsys
    #     (x::SVector{N,T}, u, tstep) ->
    #         RungeKutta4(F_sys, x, u, -tstep, nsys)::SVector{N,T}
    # end
    return VarTSystem{Nx,Nu,T}(min_tstep,max_tstep,f_sys,sys_map)
end

function NewSimpleSystemRK4(tstep::Real, f_sys::Function, nsys::Int, xtmp::SVector{Nx,T}, utmp::SVector{Nu,T}) where {Nx,Nu,T}
    sys_map = let nsys = nsys
        (x::SVector{Nx,T}, u::SVector{Nu,T}) ->
        RungeKutta4(f_sys, x, u, tstep, nsys)::SVector{Nx,T} #use Vector instead?
    end
    # sys_inv_map = let nsys = nsys
    #     (x::SVector{N,T}, u, tstep) ->
    #         RungeKutta4(F_sys, x, u, -tstep, nsys)::SVector{N,T}
    # end
    return SimpleSystem{Nx,Nu,T}(tstep,sys_map)
end

function NewSimpleDiscreteSystem(tstep::Real, sys_map::Function, xtmp::SVector{Nx,T}, utmp::SVector{Nu,T}) where {Nx,Nu,T}
    return SimpleSystem{Nx,Nu,T}(tstep,sys_map)
end
function NewSimpleDiscreteSystem(tstep::Real, sys_map::Function, xtmp::SVector{Nx,T}, utmp::T) where {Nx,Nu,T}
    return SimpleSystem{Nx,1,T}(tstep,sys_map)
end
function NewVarTDiscreteSystem(f::Function, F::Function, xtmp::SVector{Nx,T}, utmp::SVector{Nu,T}, min_tstep=0, max_tstep=Inf) where {Nx,Nu,T}
    return VarTDiscreteSystem{Nx,Nu,T}(min_tstep,max_tstep,f,F)
end

function computeTransition(sys::SimpleSystem,x,u)
    return sys.F(x,u)
end
function computeTransition(sys::AbstractSystem,x,u,tstep)
    return sys.F(x,u,tstep)
end
