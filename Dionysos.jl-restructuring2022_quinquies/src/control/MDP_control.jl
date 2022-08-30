import LinearAlgebra: norm
import Plots
import ProgressMeter

function value_iteration(r::SR.Concrete2ProbSymb_Relation,T;lambda=0.95,eps=0.01,MAXIT=1000)
    states = DO.enum_elems(r.xmap.dom2, 1) #only nodes marked as 1
    # setdiff!(states,T)

    n_states = length(states)
    println("n_states ", n_states)

    state2index = Dict(states .=> [i for i=1:n_states])
    for s in T
        state2index[s] = -1
    end
    # curr_rewards = zeros(MVector{n_states,Float64})
    curr_rewards = zeros(Float64, n_states)
    best_actions = Dict(states .=> zeros(Int,n_states))
    a = 1

    P = DO.enum_elems(r.xmap.dom2, -1) #prohibited states are labaled as -1
    for s in P
        curr_rewards[state2index[s]] = -Inf
    end

    err = Inf
    th = eps*(1-lambda)/(2lambda)
    stop = 0
    it = 0
    while stop<2
        # println("here8")
        it += 1
        print(it, " ")
        # prev_rewards = SVector(curr_rewards) #takes quite a lot of time
        prev_rewards = copy(curr_rewards)
        Threads.@threads for i = 1 : n_states
            # if i>= n_states println(i) end
            curr_rewards[i] = -Inf
            for (a,action) in r.sys2.nodes[states[i]].actions
                curr_reward = 0
                ok = true
                for (s2,transition) in action.trans
                    j = get(state2index,s2,0)
                    if j > 0
                        curr_reward += transition.prob * prev_rewards[j]
                    elseif j==-1 # target reached
                        curr_reward += transition.prob * 1
                    else # obstacle or outofbounds
                        curr_reward = -Inf
                        ok = false
                    end
                end
                curr_reward = -0.001 + lambda*curr_reward
                if curr_reward > curr_rewards[i] && ok
                    curr_rewards[i] = curr_reward
                    best_actions[states[i]] = a
                end
            end
        end
        if stop == 1
            stop = 2
        else
            err_rewards = curr_rewards-prev_rewards
            err = norm(err_rewards)
            print(round(err, digits=4), " ")
            if err < th || it>=MAX_IT-1
                stop = 1
            end
        end
    end
    println("tot it ", it)
    values = Dict(states .=> curr_rewards)
    return best_actions, values
end

function backward_induction!(r::SR.Concrete2ProbSymb_Relation,T;lambda=0.99,MAXIT=1000,var_tstep=true)
    states = DO.enum_elems(r.xmap.dom2,1) #only nodes marked as 1

    n_states = length(states)
    println("n_states ", n_states)

    state2index = Dict(states .=> [i for i=1:n_states])
    for s in T
        state2index[s] = -1
    end

    curr_costs = zeros(Float64, n_states)
    best_actions = Dict(states .=> zeros(Int,n_states))

    # x_step = 0.4 #CHANGED
    # theta_step = 3.0 #CHANGED
    # x_dot_step = 8.0/15 #CHANGED
    # theta_dot_step =  4.0 * 180/pi/15 #CHANGED
    # step = SVector(x_step,theta_step,x_dot_step,theta_dot_step) #CHANGED
    # target = SVector(x_step,theta_step,8.0,4.0 * 180/pi) #CHANGED

    FINAL_REWARD = 5e2
    EXIT_PENALTY = 5e2

    p = ProgressMeter.Progress(MAXIT, dt=1, desc="Computing policy...")
    for it = 1:MAXIT
        prev_costs = copy(curr_costs)
        Threads.@threads for i = 1 : n_states
        # for i = 1 : n_states
            curr_costs[i] = Inf
            if prev_costs[i] != Inf
                for (a,action) in r.sys2.nodes[states[i]].actions
                    curr_cost = 0
                    for (s2,transition) in action.trans
                        j = get(state2index,s2,0)
                        # if states[i] == 132564
                        #     println("s2 ", s2)
                        #     println("j ", j)
                        # end
                        if j > 0
                            if var_tstep
                                curr_cost += transition.prob * (action.tstep + lambda * prev_costs[j])
                            else
                                curr_cost += transition.prob * (1 + lambda * prev_costs[j])
                            end
                        elseif j==-1 # target reached
                            if var_tstep
                                curr_cost += transition.prob * (action.tstep + lambda * (-FINAL_REWARD))
                            else
                                curr_cost += transition.prob * (1 + lambda * (-FINAL_REWARD))
                            end
                        else # obstacle or outofbounds
                            if var_tstep
                                curr_cost += transition.prob * (action.tstep + lambda * EXIT_PENALTY)
                            else
                                curr_cost += transition.prob * (1 + lambda * EXIT_PENALTY)
                            end
                        end
                    end
                    # if curr_cost == Inf
                    #     delete!(r.sys2.nodes[states[i]].actions,a)
                    # else
                    if curr_cost < curr_costs[i]
                        curr_costs[i] = curr_cost
                        best_actions[states[i]] = a
                    end
                end
            end
            # if curr_costs[i] == Inf
            #     best_actions[states[i]] = 0
            # end
        end
        # if true
        #     costs = Dict(states .=> curr_costs) #CHANGED
        #     fig = Plots.plot(aspect_ratio = x_step/theta_step,legend = false) #CHANGED
        #     plot_costs!(r, costs,fig=fig) #CHANGED
        # end
        ProgressMeter.next!(p)
    end
    costs = Dict(states .=> curr_costs)
    return best_actions, costs
end

function compute_reachable_set(r::SR.Concrete2ProbSymb_Relation,T,mark_unreachable=false)
    visited = Set(T)
    to_visit = Set(T)
    while !isempty(to_visit)
        s2 = pop!(to_visit)
        push!(visited,s2)
        entering_nodes = Set(ST.list_entering_nodes(r.sys2,s2))
        new_nodes = setdiff(entering_nodes,visited)
        to_visit = union(to_visit,new_nodes)
    end
    if mark_unreachable
        states = Set(DO.enum_elems(r.xmap.dom2, 1))
        unreachable = setdiff(states,visited)
        for s1 in unreachable
            DO.set_s!(r.xmap.dom2,s1,-2)
        end
    end
    return visited
end

function simulate(r::SR.Concrete2ProbSymb_Relation{Nx1,Nu1,T1,Nx2,Nu2,T2},x,T,best_actions,values,fig;var_tstep=true,tstep=0,render_each_step=true) where {Nx1,Nu1,T1,Nx2,Nu2,T2}
    steps_count = 0
    ttot = 0
    s0 = 0
    s1 = 0
    s2 = 0
    s3 = 0
    s4 = 0
    u = nothing
    try
        while !in(x,T)
            steps_count += 1
            # print(x)
            s = DM.map1(r.xmap,SVector{Nx1}(x))
            s4 = s3
            s3 = s2
            s2 = s1
            s1 = s0
            s0 = s
            if s0==s3 || s0==s4
                println("encountered loop")
                println(u)
                #keeps old u
            else
                a = get(best_actions,s,0)
                u = DM.map2(r.umap,a)
            end
            # println(x)
            xprev = x
            if var_tstep
                if tstep == 0
                    tstep = r.sys2.nodes[s].actions[a].tstep
                end
                ttot += tstep
                x = ST.computeTransition(r.sys1,xprev,u,tstep)
                # println("tstep ", tstep)
                # println("h,f ", DO.get_cell_center_from_coord(r.xmap.dom1, xprev), " ", r.sys1.f(x,u))
                # h = UT.size(DO.get_smallest_cell_from_coord(r.xmap.dom1,x))
                # if (min(abs.(x.-xprev)...) < min(h./2 ...))
                #     ttot += tstep
                #     x = ST.computeTransition(r.sys1,x,u,tstep)
                # end
            else
                tstep = r.sys1.tstep
                ttot += tstep
                x = ST.computeTransition(r.sys1,xprev,u)
                # h = UT.size(DO.get_smallest_cell_from_coord(r.xmap.dom1,x))
                # if (min(abs.(x.-xprev)...) < min(h./2 ...))
                #     ttot += tstep
                #     x = ST.computeTransition(r.sys1,x,u)
                # end
            end
            # println(ttot)
            x = DO.bring_into_period(r.xmap.dom1, x) # if periodic
            # actions = ST.list_actions(r.sys2,s)
            # for a in actions
            #     println(a)
            #     println(ST.list_trans(r.sys2,s,a))
            # end
            # while !DO.isin(r.xmap.dom1,x)
            #     a = pop!(actions)
            #     u = DM.map2(r.umap,a)
            #     x = ST.computeTransition(r.sys1,xprev,u)
            # end
            Plots.plot!([[xprev[1], x[1]]], [[xprev[2], x[2]]], opacity=1, linewidth=2, color=:blue)
            if render_each_step
                display(fig)
            end
            println(x)
            # if s0==1349
            #     break
            # end
        end
        println("reached in ", steps_count, " steps")
        println("reached in ", ttot, " seconds")
    finally
        display(fig)
    end
end

function simulateRK4(RKsys,r::SR.Concrete2ProbSymb_Relation{Nx1,Nu1,T1,Nx2,Nu2,T2},x,T,best_actions,values,fig;var_tstep=true,tstep=0,render_each_step=true) where {Nx1,Nu1,T1,Nx2,Nu2,T2}
    steps_count = 0
    ttot = 0
    try
        while !in(x,T)
            steps_count += 1
            # print(x)
            s = DM.map1(r.xmap,SVector{Nx1}(x))
            # println(values[s])
            # println(s)
            # if s==27
            #     println(best_actions[27])
            # end
            a = get(best_actions,s,0)
            u = DM.map2(r.umap,a)
            # println(x)
            xprev = x
            if var_tstep
                if tstep == 0
                    tstep = r.sys2.nodes[s].actions[a].tstep
                end
                ttot += tstep
                x = ST.computeTransition(RKsys,xprev,u,tstep)
                # println("tstep ", tstep)
                # println("h,f ", DO.get_cell_center_from_coord(r.xmap.dom1, xprev), " ", r.sys1.f(x,u))
                # h = UT.size(DO.get_smallest_cell_from_coord(r.xmap.dom1,x))
                # if (min(abs.(x.-xprev)...) < min(h./2 ...))
                #     ttot += tstep
                #     x = ST.computeTransition(r.sys1,x,u,tstep)
                # end
            else
                tstep = r.sys1.tstep
                ttot += tstep
                x = ST.computeTransition(RKsys,xprev,u)
                # h = UT.size(DO.get_smallest_cell_from_coord(r.xmap.dom1,x))
                # if (min(abs.(x.-xprev)...) < min(h./2 ...))
                #     ttot += tstep
                #     x = ST.computeTransition(r.sys1,x,u)
                # end
            end
            println(ttot)
            x = DO.bring_into_period(r.xmap.dom1, x) # if periodic
            # actions = ST.list_actions(r.sys2,s)
            # for a in actions
            #     println(a)
            #     println(ST.list_trans(r.sys2,s,a))
            # end
            # while !DO.isin(r.xmap.dom1,x)
            #     a = pop!(actions)
            #     u = DM.map2(r.umap,a)
            #     x = ST.computeTransition(r.sys1,xprev,u)
            # end
            Plots.plot!([[xprev[1], x[1]]], [[xprev[2], x[2]]], opacity=1, linewidth=2, color=:blue)
            if render_each_step
                display(fig)
            end
        end
        println("reached in ", steps_count, " steps")
        println("reached in ", ttot, " seconds")
    finally
        display(fig)
    end
end

function simulate_cartpole(r::SR.Concrete2ProbSymb_Relation{Nx1,Nu1,T1,Nx2,Nu2,T2},xi,T,best_actions,values,fig) where {Nx1,Nu1,T1,Nx2,Nu2,T2}
    x = xi
    steps_count = 0
    target_reached = false
    while target_reached<1
        steps_count += 1
        println("x=",x[1]," th=",x[2])
        println("x_dot=",x[3]," th_dot=",x[4])
        s = DM.map1(r.xmap,SVector{Nx1}(x))
        # println(values[s])
        a = get(best_actions,s,0)
        if a==0
            println("notfound")
            u = SVector(0.0)
        else
            u = DM.map2(r.umap,a)
            println("value= ", values[s])
        end
        if in(x,T)
            target_reached = true
        end
        println("u=",u)
        xprev = x
        x = ST.computeTransition(r.sys1,xprev,u)
        # actions = ST.list_actions(r.sys2,s)
        # while !DO.isin(r.xmap.dom1,x)
        #     a = pop!(actions)
        #     u = DM.map2(r.umap,a)
        #     x = ST.computeTransition(r.sys1,xprev,u)
        # end
        Plots.plot!([[xprev[1], x[1]]], [[xprev[2], x[2]]], opacity=1, linewidth=2, color=:blue)
        display(fig)
    end
    println("reached in ", steps_count, " steps")
end

function shape2D(rect::UT.HyperRectangle)
    Plots.Shape([(rect.lb[1],rect.lb[2]),(rect.lb[1],rect.ub[2]),(rect.ub[1],rect.ub[2]),(rect.ub[1],rect.lb[2])])
end
# function shape_rect2D(rect::UT.HyperRectangle{SVector})
#     Plots.Shape([(rect.lb[1],rect.lb[2]),(rect.lb[1],rect.ub[2]),(rect.ub[1],rect.ub[2]),(rect.ub[1],rect.lb[2])])
# end
# CHECKED
function plot_the_domain!(dom::DO.AbstractNestedDomain{N,T}) where {N,T}
    dims = size(dom.table)
    pos = ones(Int,N)
    for i = 1 : dims[1]
        pos[1] = i
        for j = 1 : dims[2]
            pos[2] = j
            iscell = dom.table[pos...]
            if iscell == 1
                cell = DO.get_cell_from_pos(dom, SVector{N}(pos))
                Plots.plot!(shape2D(cell), opacity=0.0,color=:green)
            elseif iscell == 0
                cell = DO.get_cell_from_pos(dom, SVector{N}(pos))
                Plots.plot!(shape2D(cell), opacity=1.0,color=:red)
            else
                subdom = dom.subdomains[pos][]
                plot_the_domain!(subdom)
            end
        end
    end
end
function plot_the_values!(r,values;scale=1)
    drawn_rects = Dict{UT.HyperRectangle,Float64}()
    p = ProgressMeter.Progress(length(values), dt=1, desc="Plotting...")
    for (s,v) in values
        cell = DO.get_smallest_cell_from_coord(r.xmap.dom1, DM.map2(r.xmap,s))
        rect = UT.buildRectangle(SVector(cell.lb[1],cell.lb[2]),SVector(cell.ub[1],cell.ub[2]))
        old_v = get(drawn_rects,rect,nothing)
        if old_v==nothing || v>old_v
            if v>= 0
                Plots.plot!(shape2D(rect), opacity=1,color=:white)
                Plots.plot!(shape2D(rect), opacity=v/scale,color=:green)
            else
                Plots.plot!(shape2D(rect), opacity=1,color=:white)
                Plots.plot!(shape2D(rect), opacity=-v/scale,color=:red)
            end
            drawn_rects[rect] = v
        end
        ProgressMeter.next!(p)
    end
end
#SLOW, visits all states and plots them multiple times
function plot_values!(r::SR.Concrete2ProbSymb_Relation,values;scale=1,init=nothing, target=nothing,fig=Plots.plot(aspect_ratio = 1,legend = false))
    plot_the_domain!(r.xmap.dom1)
    plot_the_values!(r,values,scale=scale)
    if init != nothing
        Plots.plot!(shape2D(init), opacity=1,color=:orange)
    end
    if target != nothing
        Plots.plot!(shape2D(target), opacity=1,color=:blue)
    end
    display(fig)
end
function plot_costs!(r::SR.Concrete2ProbSymb_Relation,costs;init=nothing, target=nothing,fig=Plots.plot(aspect_ratio = 1,legend = false))
    values = Dict{Int,Float64}()
    for (s,c) in costs
        # values[s] = 1/c
        values[s] = -c/1e3
        # values[s] = -c
    end
    plot_values!(r,values,init=init, target=target,fig=fig)
end

function plot_the_values_at_origin!(r,values;scale=1)
    p = ProgressMeter.Progress(length(values), dt=1, desc="Plotting...")
    for (s,v) in values
        center = DM.map2(r.xmap,s)
        cell = DO.get_smallest_cell_from_coord(r.xmap.dom1, center)
        rect = UT.buildRectangle(SVector(cell.lb[1],cell.lb[2]),SVector(cell.ub[1],cell.ub[2]))
        # if ((-1.0 < center[3] < 1.0) && (-1.0 < center[4] < 1.0))
        #     println(center)
        # end
        # origin = SVector(0.0,0.0,0.0,0.0)
        # origin_cell = DO.get_smallest_cell_from_coord(r.xmap.dom1, origin)
        # println("origin", UT.center(origin_cell))
        if (isapprox(center[3], 0.0, atol=1e-13) && isapprox(center[4], 0.0, atol=1e-13))
            if v>= 0
                Plots.plot!(shape2D(rect), opacity=1,color=:white)
                Plots.plot!(shape2D(rect), opacity=v/scale,color=:green)
            else
                Plots.plot!(shape2D(rect), opacity=1,color=:white)
                Plots.plot!(shape2D(rect), opacity=-v/scale,color=:red)
            end
        end
        ProgressMeter.next!(p)
    end
end
function plot_values_at_origin!(r::SR.Concrete2ProbSymb_Relation,values;scale=1,init=nothing, target=nothing,fig=Plots.plot(aspect_ratio = 1,legend = false))
    plot_the_domain!(r.xmap.dom1)
    plot_the_values_at_origin!(r,values,scale=scale)
    if init != nothing
        Plots.plot!(shape2D(init), opacity=1,color=:orange)
    end
    if target != nothing
        Plots.plot!(shape2D(target), opacity=1,color=:blue)
    end
    display(fig)
end
function plot_costs_at_origin!(r::SR.Concrete2ProbSymb_Relation,costs;init=nothing, target=nothing,fig=Plots.plot(aspect_ratio = 1,legend = false))
    values = Dict{Int,Float64}()
    for (s,c) in costs
        # values[s] = 1/c
        values[s] = -c/5e2
        # values[s] = -c
    end
    plot_values_at_origin!(r,values,init=init, target=target,fig=fig)
end

function value_iteration_safety(r::SR.Concrete2ProbSymb_Relation,target::SVector,accuracy::SVector;lambda=0.999,MAXIT=200)
    safeSet = UT.HyperRectangle(target-accuracy, target+accuracy)
    safeStates = DM.map1(r.xmap, safeSet)
    last_safe = length(safeStates)
    # unsafeStates = setdiff(DO.enum_elems(r.xmap.dom2),safeStates)
    # states = append!(safeStates,unsafeStates)
    states = union(safeStates, DO.enum_elems(r.xmap.dom2))

    n_states = length(states)
    println("n_states ", n_states)

    state2index = Dict(states .=> [i for i=1:n_states])
    curr_rewards = zeros(Float64, n_states).-10000
    prev_rewards = zeros(Float64, n_states)
    best_actions = Dict(states .=> zeros(Int,n_states))

    p = ProgressMeter.Progress(MAXIT, dt=1, desc="Computing policy...")
    for it = 1:MAXIT
        Threads.@threads for i = 1 : n_states
        # for i = 1 : n_states
            curr_rewards[i] = -Inf
            for (a,action) in r.sys2.nodes[states[i]].actions
                curr_reward = 0
                for (s2,transition) in action.trans
                    j = state2index[s2]
                    if j<=last_safe
                        curr_reward += transition.prob * (1 + lambda*prev_rewards[j])
                    else
                        curr_reward += transition.prob * (-100 + lambda*prev_rewards[j])
                    end
                end
                if curr_reward > curr_rewards[i]
                    curr_rewards[i] = curr_reward
                    best_actions[states[i]] = a
                end
            end
        end
        ProgressMeter.next!(p)
        prev_rewards = copy(curr_rewards)
    end
    rewards = Dict(states .=> curr_rewards)
    return best_actions, rewards
end

function vi_target_cart(r::SR.Concrete2ProbSymb_Relation,T::UT.HyperRectangle,B::UT.HyperRectangle;lambda=0.99,MAXIT=500)
    targets = DM.map1(r.xmap, T)
    last_target = length(targets)
    safe_states = union(targets, DM.map1(r.xmap, B))
    last_safe = length(safe_states)
    states = union(safe_states, DO.enum_elems(r.xmap.dom2))
    # setdiff!(safeStates,targets)

    n_states = length(states)
    println("n_states ", n_states)

    state2index = Dict(states .=> [i for i=1:n_states])
    curr_rewards = zeros(Float64, n_states).-100
    # curr_rewards = fill(-Inf, n_states)
    prev_rewards = zeros(Float64, n_states)
    best_actions = Dict(states .=> zeros(Int,n_states))

    p = ProgressMeter.Progress(MAXIT, dt=1, desc="Computing policy...")
    for it = 1:MAXIT
        # fill!(curr_rewards,-Inf)
        Threads.@threads for i = 1 : n_states
        # for i = 1 : n_states
            curr_rewards[i] = -100
            for (a,action) in r.sys2.nodes[states[i]].actions
                curr_reward = 0
                for (s2,transition) in action.trans
                    j = state2index[s2]
                    if j<=last_target
                    # if in(s2,targets)
                        curr_reward += transition.prob * (1 + lambda*prev_rewards[j])
                    elseif j<=last_safe
                    # elseif in(s2,safe_states)
                        curr_reward += transition.prob * (0 + lambda*prev_rewards[j])
                    else
                        curr_reward += transition.prob * (-1 + lambda*prev_rewards[j])
                    end
                end
                if curr_reward > curr_rewards[i]
                    curr_rewards[i] = curr_reward
                    best_actions[states[i]] = a
                end
            end
        end
        ProgressMeter.next!(p)
        prev_rewards = copy(curr_rewards)
    end
    rewards = Dict(states .=> curr_rewards)
    return best_actions, rewards
end

function vi_target_cart_bis(r::SR.Concrete2ProbSymb_Relation,T::UT.HyperRectangle;MAXIT=500)
    targets = DM.map1(r.xmap, T)
    last_target = length(targets)
    states = union(targets, DO.enum_elems(r.xmap.dom2))
    # setdiff!(safeStates,targets)

    n_states = length(states)
    println("n_states ", n_states)

    state2index = Dict(states .=> [i for i=1:n_states])
    curr_rewards = zeros(Float64, n_states)
    # curr_rewards = fill(-Inf, n_states)
    prev_rewards = zeros(Float64, n_states)
    best_actions = Dict(states .=> zeros(Int,n_states))

    p = ProgressMeter.Progress(MAXIT, dt=1, desc="Computing policy...")
    for it = 1:MAXIT
        # fill!(curr_rewards,-Inf)
        Threads.@threads for i = 1 : n_states
        # for i = 1 : n_states
            curr_rewards[i] = -1
            for (a,action) in r.sys2.nodes[states[i]].actions
                curr_reward = 0
                for (s2,transition) in action.trans
                    j = state2index[s2]
                    if j<=last_target
                    # if in(s2,targets)
                        curr_reward += transition.prob * (1 + prev_rewards[j])
                    else
                        curr_reward += transition.prob * (0 + prev_rewards[j])
                    end
                end
                if curr_reward > curr_rewards[i]
                    curr_rewards[i] = curr_reward
                    best_actions[states[i]] = a
                end
            end
        end
        ProgressMeter.next!(p)
        prev_rewards = copy(curr_rewards)
    end
    rewards = Dict(states .=> curr_rewards)
    return best_actions, rewards
end

function bi_target(r::SR.Concrete2ProbSymb_Relation,T::UT.HyperRectangle;MAXIT=500)

    targets = DM.map1(r.xmap, T)
    last_target = length(targets)
    reached = copy(targets)
    just_reached = copy(targets)

    states = union(targets, DO.enum_elems(r.xmap.dom2))
    n_states = length(states)
    println("n_states ", n_states)

    state2index = Dict(states .=> [i for i=1:n_states])
    curr_rewards = zeros(Float64, n_states) .-1
    prev_rewards = zeros(Float64, n_states)
    best_actions = Dict(states .=> zeros(Int,n_states))
    is_expanding = true

    p = ProgressMeter.Progress(MAXIT, dt=1, desc="Computing policy...")
    for it = 1:MAXIT
        if is_expanding
            new_nodes = Int[]
            for s in just_reached
                new_nodes = union(new_nodes,ST.list_entering_nodes(r.sys2,s))
            end
            just_reached = setdiff(new_nodes,reached)
            reached = append!(reached,just_reached)
            n_just_reached = length(just_reached)
            if just_reached==0
                is_expanding = false
            end
        end
        curr_rewards .= -1
        Threads.@threads for s1 in reached
        # for i = 1 : n_states
            for (a,action) in r.sys2.nodes[s1].actions
                curr_reward = 0
                for (s2,transition) in action.trans
                    j = state2index[s2]
                    if j<=last_target
                    # if in(s2,targets)
                        curr_reward += transition.prob * (1 + prev_rewards[j])
                    else
                        curr_reward += transition.prob * prev_rewards[j]
                    end
                end
                if curr_reward > curr_rewards[state2index[s1]]
                    curr_rewards[state2index[s1]] = curr_reward
                    best_actions[s1] = a
                end
            end
        end
        ProgressMeter.next!(p)
        prev_rewards = copy(curr_rewards)
    end
    rewards = Dict(states .=> curr_rewards)
    return best_actions, rewards
end

function fast_bi_target(r::SR.Concrete2ProbSymb_Relation,T::UT.HyperRectangle)
    to_visit = DM.map1(r.xmap, T)
    to_visit_next = Int[]
    reached = copy(to_visit)
    states = DO.enum_elems(r.xmap.dom2)
    n_states = length(states)
    println("n_states ", n_states)

    costs = Dict(states .=> zeros(Int,n_states))
    best_actions = Dict(states .=> zeros(Int,n_states))

    visit = true
    cost = 1
    while !isempty(to_visit)
        s2 = pop!(to_visit)
        entering_t = r.sys2.nodes[s2].entering_trans
        for (s1,A) in entering_t
            if !in(s1,reached)
                costs[s1] = cost
                best_actions[s1] = maximum(A)
                push!(reached,s1)
                push!(to_visit_next,s1)
            end
        end
        if isempty(to_visit)
            append!(to_visit,to_visit_next)
            empty!(to_visit_next)
            cost += 1
        end
    end
    return best_actions, costs
end

function simulateDC(r::SR.Concrete2ProbSymb_Relation,xi,best_actions)
    x = xi
    count = 0
    while true
        count += 1
        # if count%10 == 0
        println("x ", x)
        s = DM.map1(r.xmap,x)
        println("s ", s)
        a = best_actions[s]
        u = DM.map2(r.umap,a)
        x = ST.computeTransition(r.sys1,x,u)
    end
end
