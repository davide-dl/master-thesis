import ProgressMeter

# TODO: #change abstract fields? (probably not necessary)
struct Concrete2ProbSymb_Relation{Nx1,Nu1,T1,Nx2,Nu2,T2<:Int} <: AbstractRelation{Nx1,Nu1,T1,Nx2,Nu2,T2}
    xmap::Map.AbstractDomainMap{Nx1,T1,Nx2,T2}
    umap::Map.AbstractDomainMap{Nu1,T1,Nu2,T2}
    sys1::System.AbstractSystem{Nx1,Nu1,T1}
    sys2::System.AbstractSystem{Nx2,Nu2,T2}
    samples_per_dim::Int
end

function buildConcrete2ProbSymb_Relation(xmap::Map.AbstractDomainMap, umap::Map.AbstractDomainMap, sys1::System.AbstractSystem{Nx1,Nu1,T1},sys2::System.AbstractSystem{Nx2,Nu2,T2},n::Int) where {Nx1,Nu1,T1,Nx2,Nu2,T2}
    return Concrete2ProbSymb_Relation{Nx1,Nu1,T1,Nx2,Nu2,T2}(xmap,umap,sys1,sys2,n)
end

# it's assumed that s1 maps to a HyperRectangle {x1}
# all the transitions (s1,a)->s2 are computed
# it's assumed that a maps to a single input u
# the transitions F(x1,u)->x2 are assumed to be deterministic
function computeTransitions!(r::Concrete2ProbSymb_Relation{Nx1,Nu1,T1,Nx2,Nu2,T2}, s::Int, for_all_reachable, computedNodes) where {Nx1,Nu1,T1,Nx2,Nu2,T2}
    # suppose the state is already in the domain
    tocompute = Set(s)
    n_samples = r.samples_per_dim^Nx1

    while !isempty(tocompute)
        s1 = pop!(tocompute)
        r.sys2.nodes[s1] = get(r.sys2.nodes, s1, System.buildEmptyNode())
        computedNodes[s1] = true

        center = Map.map2(r.xmap,s1)
        cell = DO.get_smallest_cell_from_coord(r.xmap.dom1,center)
        step = UT.size(cell)./(r.samples_per_dim+1)
        lb = cell.lb + step./2

        x1 = Vector{T1}(undef,Nx1)
        is_node_prohibited = true
        for a in DO.enum_elems(r.umap.dom2) #all actions are considered
            # println(a)
            u = Map.map2(r.umap,a)
            is_action_allowed = true
            action = System.buildEmptyAction()
            tmp_dict = Dict{Int64,Int64}() # state,count
            for k = 0 : n_samples - 1
                tmp = k
                for i = 1 : Nx1
                    x1[i] = (tmp % r.samples_per_dim)*step[i] + lb[i]
                    tmp = trunc(Int, tmp / r.samples_per_dim) + 1
                end
                x2_tmp = System.computeTransition(r.sys1, SVector{Nx1,T1}(x1), u)
                x2 = DO.bring_into_period(r.xmap.dom1, x2_tmp) # if periodic

                if DO.isin(r.xmap.dom1,x2)
                    s2 = Map.map1(r.xmap,x2)
                    tmp_dict[s2] = get(tmp_dict, s2, 0) + 1
                    if for_all_reachable
                        alreadyComputed = get(computedNodes, s2, false)
                        if !alreadyComputed
                            push!(tocompute,s2)
                        end
                    end
                else # node goes into an obstacle or outofbounds
                    is_action_allowed = false
                end
            end
            if is_action_allowed
                r.sys2.nodes[s1].actions[a] = action
                is_node_prohibited = false
                for (s2,count) in tmp_dict
                    System.addEnteringTrans!(r.sys2,s1,a,s2)
                    action.trans[s2] = System.buildTransition(count/n_samples)
                end
            # else
            #     println("action not allowed! a = ",a," s1 = ", s1)
            end
        end
        if is_node_prohibited
            DO.set_s!(r.xmap.dom2,s2,-1)
        end
    end
end

# TODO: use a list of rects as init
function buildFromInitSet!(r::Concrete2ProbSymb_Relation, I::UT.HyperRectangle)
    Xdom = r.xmap.dom1
    # DO.addSet!(Xdom,I)
    X = DO.enum_centers_in_rect(Xdom, I)
    computedNodes = Dict{Int,Bool}()
    for x1 in X
        s1 = Map.map1(r.xmap, x1)
        alreadyComputed = get(computedNodes, s1, false)
        if !alreadyComputed
            computeTransitions!(r,s1,true, computedNodes)
        end
    end
    allStates = DO.enum_elems(r.xmap.dom2)
    for s in allStates
        alreadyComputed = get(computedNodes, s, false)
        if !alreadyComputed
            DO.set_s!(r.xmap.dom2,s,-2) #unreachable
        end
    end
end

function compute_avg_tstep(r, center, h, step, lb, n_samples, x1, u, Nx1, T1) #x1 is empty
    sum_tstep = 0
    for k = 0 : n_samples - 1
        tmp = k
        for i = 1 : Nx1
            x1[i] = (tmp % r.samples_per_dim)*step[i] + lb[i]
            tmp = trunc(Int, tmp / r.samples_per_dim) + 1
        end
        x1_dot = r.sys1.f(SVector{Nx1,T1}(x1), u)
        d = h - sign.(x1_dot).*(center.-x1)
        sum_tstep += min(abs.(d./x1_dot)...)
    end
    tstep = sum_tstep/n_samples
    tstep = max(tstep,r.sys1.min_tstep)
    return min(tstep,r.sys1.max_tstep)
end
function compute_central_tstep(r, center, h, u) #x1 is empty
    x_dot = r.sys1.f(center, u)
    tstep = min(abs.(h./x_dot)...)
    tstep = max(tstep,r.sys1.min_tstep)
    return min(tstep,r.sys1.max_tstep)
end
function compute_central_tstep_neighbor(r, center, h, u, Nx1, T1) #x1 is empty
    x_dot = r.sys1.f(center, u)
    direction = sign.(x_dot)
    h_nb = Vector{T1}(undef,Nx1) #h neighbors in that direction
    for i = 1 : Nx1
        x_nb = Vector(center)
        x_nb[i] += direction[i]*(2*eps(T1)+h[i]/2) #just outside the border on the given direction
        isin, nb_dom, nb_pos = DO.get_smallest_pos_from_coord(r.xmap.dom1,SVector{Nx1,T1}(x_nb))
        if isin
            cell_nb = DO.get_cell_from_pos(nb_dom,nb_pos)
            h_nb[i] = UT.size(cell_nb)[i]
        else
            h_nb[i] = Inf
        end
    end
    tstep = min(abs.( (h.+h_nb)./2 ./x_dot)...)
    tstep = max(tstep,r.sys1.min_tstep)
    return min(tstep,r.sys1.max_tstep)
end
function compute_avg_tstep_neighbor(r, center, h, step, lb, n_samples, x1, u, Nx1, T1) #x1 is empty
    h_nb = Vector{T1}(undef,Nx1) #h neighbors in that direction
    sum_tstep = 0
    for k = 0 : n_samples - 1

        #compute x
        tmp = k
        for i = 1 : Nx1
            x1[i] = (tmp % r.samples_per_dim)*step[i] + lb[i]
            tmp = trunc(Int, tmp / r.samples_per_dim) + 1
        end
        x_dot = r.sys1.f(SVector{Nx1,T1}(x1), u)

        direction = sign.(x_dot)
        for i = 1 : Nx1
            x_nb = Vector(center)
            x_nb[i] += direction[i]*(2*eps(T1)+h[i]/2) #just outside the border on the given direction
            isin, nb_dom, nb_pos = DO.get_smallest_pos_from_coord(r.xmap.dom1,SVector{Nx1,T1}(x_nb))
            if isin
                cell_nb = DO.get_cell_from_pos(nb_dom,nb_pos)
                h_nb[i] = UT.size(cell_nb)[i]
            else
                h_nb[i] = Inf
            end
        end
        d = (h.+h_nb)./2 - sign.(x_dot).*(center.-x1)
        sum_tstep = sum_tstep + min(abs.(d./x_dot)...)
    end
    tstep = sum_tstep/n_samples
    tstep = max(tstep,r.sys1.min_tstep)
    return min(tstep,r.sys1.max_tstep)
end

function computeSingleTransitions!(r::Concrete2ProbSymb_Relation{Nx1,Nu1,T1,Nx2,Nu2,T2}, s1::Int; simulate_borders=false, look_at_corners=false, var_tstep=1) where {Nx1,Nu1,T1,Nx2,Nu2,T2}
    # suppose the state is already in the domain
    # if n_samples == 0
    n_samples = r.samples_per_dim^Nx1
    # end

    # r.sys2.nodes[s1] = System.buildEmptyNode()

    center = Map.map2(r.xmap,s1)
    cell = DO.get_smallest_cell_from_coord(r.xmap.dom1,center)
    h = UT.size(cell)
    # println(cell)
    # println(center)
    # sleep(1.0)
    if simulate_borders
        step = h./(r.samples_per_dim-1)
        lb = cell.lb
    else
        step = h./(r.samples_per_dim)
        lb = cell.lb + step./2
    end

    if look_at_corners
        corners = UT.get_corners(cell)
    end

    x1 = Vector{T1}(undef,Nx1)
    is_node_prohibited = true
    #TODO: following line can be computed outside
    for a in DO.enum_elems(r.umap.dom2) #all actions are considered
        u = Map.map2(r.umap,a)
        is_action_allowed = true
        tmp_dict = Dict{Int64,Int64}() # state,count

        if var_tstep > 0
            if var_tstep == 1
                tstep = compute_central_tstep(r, center, h, u)
            elseif var_tstep == 2
                tstep = compute_avg_tstep(r, center, h, step, lb, n_samples, x1, u, Nx1, T1)
            elseif var_tstep == 3
                tstep = compute_central_tstep_neighbor(r, center, h, u, Nx1, T1)
            else
                tstep = compute_avg_tstep_neighbor(r, center, h, step, lb, n_samples, x1, u, Nx1, T1)
            end
            action = System.buildEmptyActionT(tstep)
        else
            action = System.buildEmptyAction()
        end


        for k = 0 : n_samples - 1
            tmp = k
            for i = 1 : Nx1
                x1[i] = (tmp % r.samples_per_dim)*step[i] + lb[i]
                tmp = trunc(Int, tmp / r.samples_per_dim) + 1
            end
            if var_tstep>0
                x2_tmp = System.computeTransition(r.sys1, SVector{Nx1,T1}(x1), u, tstep)
            else
                x2_tmp = System.computeTransition(r.sys1, SVector{Nx1,T1}(x1), u)
            end
            x2 = DO.bring_into_period(r.xmap.dom1, x2_tmp) # if periodic
            if DO.isin(r.xmap.dom1,x2)
                s2 = Map.map1(r.xmap,x2)
                tmp_dict[s2] = get(tmp_dict, s2, 0) + 1
            else # node goes into an obstacle or outofbounds
                is_action_allowed = false
                break
            end
        end
        if look_at_corners && is_action_allowed
            for x1 in corners
                if var_tstep>0
                    x2_tmp = System.computeTransition(r.sys1, SVector{Nx1,T1}(x1), u, r.sys1.min_tstep)
                    x2 = DO.bring_into_period(r.xmap.dom1, x2_tmp) # if periodic
                    if !DO.isin(r.xmap.dom1,x2)
                        is_action_allowed = false
                    end
                    x2_tmp = System.computeTransition(r.sys1, SVector{Nx1,T1}(x1), u, tstep)
                else
                    x2_tmp = System.computeTransition(r.sys1, SVector{Nx1,T1}(x1), u)
                end
                x2 = DO.bring_into_period(r.xmap.dom1, x2_tmp) # if periodic
                if !DO.isin(r.xmap.dom1,x2)
                    is_action_allowed = false
                end
            end
        end
        if is_action_allowed
            node = get(r.sys2.nodes,s1,nothing)
            while node==nothing
                sleep(0.01)
                node = get(r.sys2.nodes,s1,nothing)
                print("sleep")
            end
            is_node_prohibited = false
            r.sys2.nodes[s1].actions[a] = action
            for (s2,count) in tmp_dict
                System.addEnteringTrans!(r.sys2,s1,a,s2)
                action.trans[s2] = System.buildTransition(count/n_samples)
            end
        else
            delete!(r.sys2.nodes[s1].actions,a)
        end
    end
    if is_node_prohibited
        DO.set_s!(r.xmap.dom2,s1,-1)
        # println("prohibited", Map.map2(r.xmap,s1))
    end
    # if s1==2887
    #     println("is_node_prohibited ", is_node_prohibited)
    # end
end
function buildAll!(r::Concrete2ProbSymb_Relation{Nx1,Nu1,T1,Nx2,Nu2,T2}; simulate_borders=false, look_at_corners=false, var_tstep=1) where {Nx1,Nu1,T1,Nx2,Nu2,T2}
    S = DO.enum_elems(r.xmap.dom2)
    p = ProgressMeter.Progress(length(S), dt=1, desc="Initializing abstraction...")
    # Threads.@threads for s1 in S
    for s1 in S
        r.sys2.nodes[s1] = System.buildEmptyNode()
        ProgressMeter.next!(p)
    end
    p = ProgressMeter.Progress(length(S), dt=1, desc="Building abstraction...")
    Threads.@threads for s1 in S
    # for s1 in S
        # n_samples = r.samples_per_dim^Nx1
        computeSingleTransitions!(r,s1,simulate_borders=simulate_borders,look_at_corners=look_at_corners,var_tstep=var_tstep)
        ProgressMeter.next!(p)
    end
end


# function computeTrans_easy!(r,s1)
#     x1 = Map.map2(r.xmap,s1)
#     is_node_prohibited = true
#     #TODO: following line can be computed outside
#     for a in DO.enum_elems(r.umap.dom2) #all actions are considered
#         u = Map.map2(r.umap,a)
#         is_action_allowed = true
#         action = System.buildEmptyAction()
#         x2 = System.computeTransition(r.sys1, x1, u)
#         if DO.isin(r.xmap.dom1,x2)
#             s2 = Map.map1(r.xmap,x2)
#         else # node goes into an obstacle or outofbounds
#             is_action_allowed = false
#             break
#         end
#         if is_action_allowed
#             node = r.sys2.nodes[s1]
#             is_node_prohibited = false
#             r.sys2.nodes[s1].actions[a] = action
#             System.addEnteringTrans!(r.sys2,s1,a,s2)
#             action.trans[s2] = System.buildTransition(1.0)
#         end
#     end
#     if is_node_prohibited
#         DO.set_s!(r.xmap.dom2,s1,-1)
#     end
# end
# function buildAll_easy!(r::Concrete2ProbSymb_Relation{Nx1,Nu1,T1,Nx2,Nu2,T2}) where {Nx1,Nu1,T1,Nx2,Nu2,T2}
#     S = DO.enum_elems(r.xmap.dom2)
#     p = ProgressMeter.Progress(length(S), dt=1, desc="Initializing abstraction...")
#     # println("here")
#     # Threads.@threads for s1 in S
#     for s1 in S
#         r.sys2.nodes[s1] = System.buildEmptyNode()
#         ProgressMeter.next!(p)
#     end
#     p = ProgressMeter.Progress(length(S), dt=1, desc="Building abstraction...")
#     Threads.@threads for s1 in S
#     # for s1 in S
#         # n_samples = r.samples_per_dim^Nx1
#         computeTrans_easy!(r,s1)
#         ProgressMeter.next!(p)
#     end
# end
