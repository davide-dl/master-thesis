# TODO:

# UT.leaq(a,b) = (a <= b) || isapprox(a,b,atol=1e-8)
# UT.greaq(a,b) = (a >= b) || isapprox(a,b,atol=1e-8)

import Plots

"""
    AbstractNestedDomain{N,T<:Real} <: DomainType{N,T}

Description
"""
abstract type AbstractNestedDomain{N,T<:Real} <: DomainType{N,T} end # It could be a "normal" nested domain or a ObstacleNestedDomain or a FreeNestedDomain

"""
    NestedDomain{N,T} <: AbstractNestedDomain{N,T}

Description
"""
struct NestedDomain{N,T} <: AbstractNestedDomain{N,T}
    limits::UT.HyperRectangle # limits of the domain
    h::SVector{N,T} # step-size
    depth::Int #depth of the domain in the overall nested domain "tree"
                # =0 if isroot
    isroot::Bool #if it's the first big domain
    isleaf::Bool #if all its elements are actual cells
    isperiodic::Bool
    periods::SVector{N,SVector{2,T}} # start and length of the period for each dimension, [0,0] if not periodic
    updomain::Ref # reference to the "upper" domain
    updomain_pos::SVector{N,Int} #position of this domain in the updomain
    #check # obstacles::RectSet #actual obstacles, useful especially if isroot (so it's possible to remove them)
    table::Array{Int,N} # table that for each cell has a 1 if the cell is there, a 0 if the cell is not in the domain
                      # a -1 if the cell is split
    subdomains::Dict{SVector{N,Int}, Ref} # each element of the table that present a -1 is associated to a subdomain
                                                # empty if isroot
    # cells::Dict{SVector{N,Int}, Int} # each element of the table that present a 1 is associated to a symbol in the automaton
end

"""
    VariableStepDomain{N,T} <: AbstractNestedDomain{N,T}

Description
"""
struct VariableStepDomain{N,T} <: AbstractNestedDomain{N,T}
    limits::UT.HyperRectangle # limits of the domain, can be unbounded
    is_variable::SVector{N,Bool} #is the dimension of variable step_size?
    h::SVector{N,T} # step-size (if isvariable =0)
    d::Vector{SVector} #position of the cells (first_element=limits.lb, last_element=limits.ub) if isvariable, else = [0]
    depth::Int #depth of the domain in the overall nested domain "tree"
                # =0 if isroot
    isroot::Bool #if it's the first big domain
    isleaf::Bool #if all its elements are actual cells
    isperiodic::Bool
    periods::SVector{N,SVector{2,T}} # start and length of the period for each dimension, [0,0] if not periodic
    updomain::Ref # reference to the "upper" domain
    updomain_pos::SVector{N,Int} #position of this domain in the updomain
    #check # obstacles::RectSet #actual obstacles, useful especially if isroot (so it's possible to remove them)
    table::Array{Int,N} # table that for each cell has a 1 if the cell is there, a 0 if the cell is not in the domain
                      # a -1 if the cell is split
    subdomains::Dict{SVector{N,Int}, Ref} # each element of the table that present a -1 is associated to a subdomain
                                                # empty if isroot
    # cells::Dict{SVector{N,Int}, Int} # each element of the table that present a 1 is associated to a symbol in the automaton
end

# this struct is useful for example if I want to remove an obstacle that intersects a cell only partially
# suppose the cell is divided in 4 eg:  |a|a|b|b|
#                                       |c|c|o|o|
#                                       |c|c|o|o|
# if "O" is the obstacle the cell is split in cells A,B,C,O
# if we are in 2D it splits the cell into 4 (or 9 if the obstacle is fully contained in the cell)

"""
    ObstacleNestedDomain{N,T} <: AbstractNestedDomain{N,T}

Description
"""
struct ObstacleNestedDomain{N,T} <: AbstractNestedDomain{N,T}
    limits::UT.HyperRectangle # limits of the domain
    depth::Int #depth of the domain in the overall nested domain "tree"
    isroot::Bool #if it's the first big domain
    isleaf::Bool #if all its elements are actual cells
    isperiodic::Bool
    periods::SVector{N,SVector{2,T}} # start and length of the period for each dimension, [0,0] if not periodic
    updomain::Ref # reference to the "upper" domain
    updomain_pos::SVector{N,Int} #position of this domain in the updomain
    obstacle::UT.HyperRectangle #it is assumed obstacle included in limits
    obstacle_pos::SVector{N,Int}
    table::Array{Int,N} # the table has 1, 2 or 3 elements for each dimension
    subdomains::Dict{SVector{N,Int}, Ref} # each element of the table that present a -1 is associated to a subdomain
    # cells::Dict{SVector{N,Int}, Int} # each element of the table that presents a 1 is associated to a symbol in the automaton
end
# this struct is useful to split a cell wrt a point
# it is useful for example when I want to delimitate precisely a certain area (eg. init set)

"""
    FreeNestedDomain{N,T} <: AbstractNestedDomain{N,T}

Description
"""
struct FreeNestedDomain{N,T} <: AbstractNestedDomain{N,T}
    limits::UT.HyperRectangle # limits of the domain
    depth::Int #depth of the domain in the overall nested domain "tree"
    isroot::Bool #if it's the first big domain
    isleaf::Bool #if all its elements are actual cells
    isperiodic::Bool
    periods::SVector{N,SVector{2,T}} # start and length of the period for each dimension, [0,0] if not periodic
    updomain::Ref # reference to the "upper" domain
    updomain_pos::SVector{N,Int} #position of this domain in the updomain
    center::SVector{N,T}
    table::Array{Int,N} # the table has 1, 2 or 3 elements for each dimension
    subdomains::Dict{SVector{N,Int}, Ref} # each element of the table that present a -1 is associated to a subdomain
    # cells::Dict{SVector{N,Int}, Int} # each element of the table that presents a 1 is associated to a symbol in the automaton
end

# CHECKED
# builds a uniform periodic domain
# for each dimension the period is specified
# 0 means non-periodic
# the limits are assumed to be inside the period
# NOTE: if a dimension is periodic is still fitted even if fit=false (this could be changed)
"""
    buildPeriodicUniformDomain(limits::UT.HyperRectangle,h::SVector{N,T},periods::SVector{N,SVector{2,T}},fit=ones(SVector{N,Bool})::SVector{N,Bool})

Description
"""
function buildPeriodicUniformDomain(limits::UT.HyperRectangle,h::SVector{N,T},periods::SVector{N,SVector{2,T}},fit=ones(SVector{N,Bool})::SVector{N,Bool}) where {N,T}
    depth = 1
    isroot = true
    isleaf = true
    updomain = Ref{AbstractNestedDomain{N,T}}()
    updomain_pos = zeros(SVector{N,Int})
    isperiodic = true

    dims = Vector{Int}(undef, N) #dimensions of the table

    newh = Vector{T}(undef,N)
    newSize = Vector{T}(undef,N)
    for i = 1:N
        if fit[i] || periods[i][0] !=0
            dims[i] = round(Int, UT.size(limits)[i]/h[i])
            newh[i] = UT.size(limits)[i]/dims[i]
            newSize[i] = UT.size(limits)[i]
        else
            dims[i] = trunc(Int, UT.size(limits)[i]/h[i])
            newSize[i] = dims[i]*h[i]
            newh[i] = h[i]
        end
    end
    h = SVector{N}(newh)
    limits = UT.buildRectangle2(limits.lb, SVector{N}(newSize))

    table = ones(Int,dims...)

    subdomains = Dict{SVector{N,Int}, Ref}() # empty since all cells are labeled as 1

    # cells = Dict{SVector{N,Int}, Int}() # empty for now. It can be constructed together with the automaton.
    #                                 # it could also be computed rn but i dunno if it's the best thing to do
    return NestedDomain(limits,h,depth,isroot,isleaf,isperiodic,periods,updomain,updomain_pos,table,subdomains)
end
# CHECKED
# builds a uniform domain
# for each dimension: (fit could be changed into a vector{bool})
#   if fit, h is adjusted wrt the limits
#   if !fit, the end of the domain is ignored

"""
    buildUniformDomain(limits::UT.HyperRectangle,h::SVector{N,T},fit=ones(SVector{N,Bool}))

Description
"""
function buildUniformDomain(limits::UT.HyperRectangle,h::SVector{N,T},fit=ones(SVector{N,Bool})) where {N,T}
    depth = 1
    isroot = true
    isleaf = true
    updomain = Ref{AbstractNestedDomain{N,T}}()
    updomain_pos = zeros(SVector{N,Int})
    isperiodic = false
    periods=zeros(SVector{N,SVector{2,T}})

    dims = Vector{Int}(undef, N)

    newh = Vector{T}(undef,N)
    newSize = Vector{T}(undef,N)
    for i = 1:N
        if fit[i]
            dims[i] = round(Int, UT.size(limits)[i]/h[i])
            newh[i] = UT.size(limits)[i]/dims[i]
            newSize[i] = UT.size(limits)[i]
        else
            dims[i] = trunc(Int, UT.size(limits)[i]/h[i])
            newSize[i] = dims[i]*h[i]
            newh[i] = h[i]
        end
    end
    h = SVector{N}(newh)
    limits = UT.buildRectangle2(limits.lb, SVector{N}(newSize))
    table = ones(Int,dims...)

    subdomains = Dict{SVector{N,Int}, Ref}() # empty since all cells are labeled as 1

    # cells = Dict{SVector{N,Int}, Int}() # empty for now. It can be constructed together with the automaton.
    #                                 # it could also be computed rn but i dunno if it's the best thing to do
    return NestedDomain(limits,h,depth,isroot,isleaf,isperiodic,periods,updomain,updomain_pos,table,subdomains)
end

"""
    buildVariableStepDomain(limits::UT.HyperRectangle,is_variable::SVector{N,Bool},h::SVector{N,T},d::Vector{SVector},periods=zeros(SVector{N,SVector{2,T}}),fit=ones(SVector{N,Bool})::SVector{N,Bool})

Description
"""
function buildVariableStepDomain(limits::UT.HyperRectangle,is_variable::SVector{N,Bool},h::SVector{N,T},d::Vector{SVector},periods=zeros(SVector{N,SVector{2,T}}),fit=ones(SVector{N,Bool})::SVector{N,Bool}) where {N,T}
    depth = 1
    isroot = true
    isleaf = true
    updomain = Ref{AbstractNestedDomain{N,T}}()
    updomain_pos = zeros(SVector{N,Int})
    isperiodic = false

    dims = Vector{Int}(undef, N) #dimensions of the table

    newh = Vector{T}(undef,N)
    newub = Vector{T}(undef,N)
    for i = 1:N
        if !is_variable[i]
            if periods[i][0] !=0
                isperiodic = true
            end
            if fit[i] || periods[i][0] !=0
                dims[i] = round(Int, UT.size(limits)[i]/h[i])
                newh[i] = UT.size(limits)[i]/dims[i]
                newub[i] = limits.ub[i]
            else
                dims[i] = trunc(Int, UT.size(limits)[i]/h[i])
                newSize[i] = limits.lb[i]+dims[i]*h[i]
                newh[i] = h[i]
            end
        else
            dims[i] = length(d[i])-1
            newub[i] = limits.ub[i]
            newh[i] = 0
        end
    end
    h = SVector{N}(newh)
    limits = UT.buildRectangle(limits.lb, SVector{N}(newub))

    table = ones(Int,dims...)

    subdomains = Dict{SVector{N,Int}, Ref}() # empty since all cells are labeled as 1

    # cells = Dict{SVector{N,Int}, Int}() # empty for now. It can be constructed together with the automaton.
    #                                 # it could also be computed rn but i dunno if it's the best thing to do
    return VariableStepDomain(limits,is_variable,h,d,depth,isroot,isleaf,isperiodic,periods,updomain,updomain_pos,table,subdomains)
end

"""
    not sure if export
"""
function computeExponentialD(center::T,smallest_step_size::T,doubles_at::T; max_cells=10, left=true, right=true, last=Inf, border=Inf) where T #exponential discretization of the interval, border always positive
    tmp = doubles_at/smallest_step_size
    e = (2+tmp) / (1+tmp)
    x = smallest_step_size
    d = [0.0]
    push!(d,x)
    for i = 1:max_cells-1
        x = smallest_step_size*e^(i)
        if UT.greaq(x,last)
            break
        end
        push!(d,x)
    end
    push!(d,border)
    if right==true && left==false
        d = d .+ center
    elseif left==true && right==false
        d = -d .+ center
    else
        r = d .+ center
        d = -d .+ center
        pop!(d)
        append!(d,r)
    end
    return d
end

# CHECKED
# return the coordinate of the small cell given the position in the big cell
# return a rectangle
"""
    get_cell_from_pos(dom::AbstractNestedDomain{N,T}, position::SVector{N,Int})

Description
"""
function get_cell_from_pos(dom::AbstractNestedDomain{N,T}, position::SVector{N,Int}) where {N,T}
    return _get_cell_from_pos(dom,position)
end
function _get_cell_from_pos(dom::NestedDomain{N,T}, position::SVector{N,Int}) where {N,T}  #tocheck
    lb = dom.limits.lb + (position.-1).*dom.h
    size = dom.h
    return UT.buildRectangle2(lb,size)
end
function _get_cell_from_pos(dom::ObstacleNestedDomain{N,T}, pos::SVector{N,Int}) where {N,T} #tocheck
    cell_lb = Vector{T}(undef,N) #empty
    cell_ub = Vector{T}(undef,N) #empty
    limits_lb = dom.limits.lb
    limits_ub = dom.limits.ub
    obstacle_lb = dom.obstacle.lb
    obstacle_ub = dom.obstacle.ub
    for i = 1:N
        if pos[i] == 1
            cell_lb[i] = limits_lb[i]
            cell_ub[i] = (dom.obstacle_pos[i] == 1) ? obstacle_ub[i] : obstacle_lb[i]
        elseif pos[i] == 2
            if dom.obstacle_pos[i] == 1
                cell_lb[i] = obstacle_ub[i]
                cell_ub[i] = limits_ub[i]
            else #dom.obstacle_pos[i] == 2
                cell_lb[i] = obstacle_lb[i]
                cell_ub[i] = obstacle_ub[i]
            end
        else # pos[i] == 3
            cell_lb[i] = obstacle_ub[i]
            cell_ub[i] = limits_ub[i]
        end
    end
    return UT.HyperRectangle(SVector{N}(cell_lb),SVector{N}(cell_ub))
end
function _get_cell_from_pos(dom::FreeNestedDomain{N,T}, pos::SVector{N,Int}) where {N,T} #tocheck
    cell_lb = Vector{T}(undef,N) #empty
    cell_ub = Vector{T}(undef,N) #empty
    for i = 1:N
        if pos[i] == 1
            cell_lb[i] = dom.limits.lb[i]
            cell_ub[i] = (dom.center[i] == dom.limits.lb[i]) ? dom.limits.ub[i] : dom.center[i]
        else
            cell_lb[i] = dom.center[i]
            cell_ub[i] = dom.limits.ub[i]
        end
    end
    return UT.HyperRectangle(SVector{N}(cell_lb),SVector{N}(cell_ub))
end
function _get_cell_from_pos(dom::VariableStepDomain{N,T}, position::SVector{N,Int}) where {N,T}  #tocheck
    lb = Vector{T}(undef,N)
    ub = Vector{T}(undef,N)
    for i=1:N
        if !is_variable[i]
            lb[i] = dom.limits.lb[i] + (position[i]-1)*dom.h[i]
            ub[i] = lb[i] + dom.h[i]
        else
            lb[i] = h[position[i]]
            ub[i] = h[position[i]+1]
        end
    end
    return UT.buildRectangle(lb,ub)
end

# I THINK IT WORKS
# returns bool,dom,pos
#returns true if the coordinate (=x) is in the domain, false if not
#if true returns also the lower-level domain and the cell position
# returns isindomain, domain, position
"""
    isin(dom::AbstractNestedDomain, coord::SVector, check_if_outofbouds=true, bringinperiod=true)

Description
"""
function isin(dom::AbstractNestedDomain, coord::SVector, check_if_outofbouds=true, bringinperiod=true)
    result,_,_ = get_smallest_pos_from_coord(dom, coord, check_if_outofbouds, bringinperiod)
    return result
end

"""
    get_smallest_pos_from_coord(dom::AbstractNestedDomain, coord::SVector, check_if_outofbouds=true, bringinperiod=true; get_lower_one=false)

Description
"""
function get_smallest_pos_from_coord(dom::AbstractNestedDomain, coord::SVector, check_if_outofbouds=true, bringinperiod=true; get_lower_one=false)
    is_in_domain, pos = get_pos_from_coord(dom, coord, check_if_outofbouds,bringinperiod, get_lower_one=get_lower_one)
    # check_if_outofbouds is used to improve performances
    if is_in_domain == 0
        if check_if_outofbouds & (pos == nothing) #outofbounds
            return false, nothing, nothing
        else
            return false, dom, pos #cell not in the the domain
        end
    elseif is_in_domain == 1
        return true, dom, pos
    else #this means that is_in_domain == -1
        return get_smallest_pos_from_coord(dom.subdomains[pos][],coord,false,false, get_lower_one=get_lower_one)
    end
end
# I THINK IT WORKS
# returns 1 if x=coord is in the current domain
#         0 if the cell is not in the domain
#         -1 if the cell is in a smallest domain
"""
    get_pos_from_coord(dom::AbstractNestedDomain{N,T}, x::SVector, check_if_outofbouds=true, bringinperiod=true; get_lower_one=false)

Description
"""
function get_pos_from_coord(dom::AbstractNestedDomain{N,T}, coord::SVector, check_if_outofbouds=true, bringinperiod=true; get_lower_one=false) where {N,T}
    return _get_pos_from_coord(dom, coord, check_if_outofbouds, bringinperiod, get_lower_one=get_lower_one)
end
function _get_pos_from_coord(dom::NestedDomain{N,T}, x::SVector, check_if_outofbouds=true, bringinperiod=true; get_lower_one=false) where {N,T}
    if bringinperiod
        coord = bring_into_period(dom,  x)
    else
        coord = x
    end
    if check_if_outofbouds && !in(coord,dom.limits) #x is outofbounds (check_if_outofbouds to improve previous function)
        return 0, nothing
    end

    if get_lower_one
        pos = Vector(ceil.(Int, (coord-dom.limits.lb) ./ dom.h))
    else
        pos = Vector(trunc.(Int, (coord-dom.limits.lb) ./ dom.h) .+ 1)
    end

    for i = 1 : N #if its on the upper bound
        if pos[i] == size(dom.table)[i] + 1
            pos[i] -= 1
        end
    end

    table_value = dom.table[pos...]
    return table_value, SVector{N}(pos)
end
# CHECKED
# returns 1 if x=coord is in the current domain
#         0 if the cell is not in the domain
#         -1 if the cell is in a smallest domain
# reminder: by construction it is assumed obstacle included in limits
function _get_pos_from_coord(dom::ObstacleNestedDomain, x::SVector{N,T}, check_if_outofbouds=true, bringinperiod=true; get_lower_one=false) where {N,T}
    if bringinperiod
        coord = bring_into_period(dom, x)
    else
        coord = x
    end
    if check_if_outofbouds && !in(coord,dom.limits) #x is outofbounds (check_if_outofbouds to improve previous function)
        return 0, nothing
    end
    dims = size(dom.table)
    pos = Vector{Int}(undef, N)
    for i = 1 : N
        if dims[i] == 1
            pos[i] = 1
        elseif dims[i] == 2
            if dom.obstacle_pos[i] == 1
                mid_point = dom.obstacle.ub[i]
            else
                mid_point = dom.obstacle.lb[i]
            end
            if get_lower_one && coord[i] ≈ mid_point
                pos[i] = 1
            else
                pos[i] = UT.greaq(coord[i],mid_point) ? 2 : 1
            end
        else #dims[i] == 3
            if get_lower_one
                pos[i] = !UT.leaq(coord[i],dom.obstacle.lb[i]) ? (!UT.leaq(coord[i],dom.obstacle.ub[i]) ? 3 : 2) : 1
            else
                pos[i] = UT.greaq(coord[i],dom.obstacle.lb[i]) ? (UT.greaq(coord[i],dom.obstacle.ub[i]) ? 3 : 2) : 1
            end
        end
    end

    table_value = dom.table[pos...]
    return table_value, SVector{N}(pos)
end
function _get_pos_from_coord(dom::FreeNestedDomain, x::SVector{N,T}, check_if_outofbouds=true, bringinperiod=true; get_lower_one=false) where {N,T}
    if bringinperiod
        coord = bring_into_period(dom, x)
    else
        coord = x
    end
    if check_if_outofbouds && !in(coord,dom.limits) #x is outofbounds (check_if_outofbouds to improve previous function)
        return 0, nothing
    end

    dims = size(dom.table)

    pos = Vector{Int}(undef, N)
    for i = 1 : N #if its on the upper bound
        if get_lower_one
            if UT.leaq(x[i],dom.center[i]) || dims[i] == 1
                pos[i] = 1
            else
                pos[i] = 2
            end
        else
            if !(UT.greaq(x[i],dom.center[i])) || dims[i] == 1
                pos[i] = 1
            else
                pos[i] = 2
            end
        end
    end

    table_value = dom.table[pos...]
    return table_value, SVector{N}(pos)
end
function _find_index(e,d::SVector)
    low = 1
    high = Size(d)

    while (low < high)
        mid = (low + high) >>> 1 #binary shift, ie. /2
        if d[mid] < e
            low = mid + 1
        else
            high = mid
        end
    end
    return mid #returns lower_index
end
function _get_pos_from_coord(dom::VariableStepDomain{N,T}, x::SVector, check_if_outofbouds=true, bringinperiod=true; get_lower_one=false) where {N,T}
    if bringinperiod
        coord = bring_into_period(dom,  x)
    else
        coord = x
    end
    if check_if_outofbouds && !in(coord,dom.limits) #x is outofbounds (check_if_outofbouds to improve previous function)
        return 0, nothing
    end

    pos = Vector{Int}(undef,N)
    if get_lower_one
        for i = 1:N
            if !is_variable[i]
                pos[i] = ceil(Int, (coord[i]-dom.limits.lb[i]) / dom.h[i])
                if pos[i] == size(dom.table)[i] + 1
                    pos[i] -= 1
                end
            else
                pos[i] = _find_index(coord[i], d[i])
                if (d[i][pos[i]] == coord[i]) && (pos[i] != size(dom.table)[i]+1)
                    pos[i] -= 1
                end
            end
        end
    else
        for i = 1:N
            if !is_variable[i]
                pos[i] = trunc(Int, (coord[i]-dom.limits.lb[i]) / dom.h[i] + 1)
                if pos[i] == size(dom.table)[i] + 1
                    pos[i] -= 1
                end
            else
                pos[i] = _find_index(coord[i], d[i])
                if (d[i][pos[i]] == coord[i]) && (pos[i] != size(dom.table)[i]+1)
                    pos[i] -= 1
                end
            end
        end
    end

    table_value = dom.table[pos...]
    return table_value, SVector{N}(pos)
end

"""
    get_smallest_cell_from_coord(dom::AbstractNestedDomain, coord::SVector, check_if_outofbouds=true, bringinperiod=true; get_lower_one=false)

Description
"""
function get_smallest_cell_from_coord(dom::AbstractNestedDomain, coord::SVector, check_if_outofbouds=true, bringinperiod=true; get_lower_one=false)
    _, subdom, pos = get_smallest_pos_from_coord(dom, coord, check_if_outofbouds, bringinperiod, get_lower_one=get_lower_one)
    cell = get_cell_from_pos(subdom,pos)
    return cell
end

# looks in what cells x (=coord) is and returns the center of the cell
"""
    get_cell_center_from_coord(dom::AbstractNestedDomain, coord::SVector, check_if_outofbouds=true, bringinperiod=true; get_lower_one=false)

Description
"""
function get_cell_center_from_coord(dom::AbstractNestedDomain, coord::SVector, check_if_outofbouds=true, bringinperiod=true; get_lower_one=false)
    cell = get_smallest_cell_from_coord(dom, coord, check_if_outofbouds, bringinperiod, get_lower_one=get_lower_one)
    return UT.center(cell)
end


# CHECKED
"""
    removeCoordFromDomain!(dom::AbstractNestedDomain,pos::SVector)

Description
"""
function removeCoordFromDomain!(dom::AbstractNestedDomain,pos::SVector)
    dom.table[pos...] = 0
    return
end

# I THINK IT WORKS
# TODO: decide what to do if isin==-1
"""
    refineCell!(dom::AbstractNestedDomain{N,T},uppos::SVector{N,Int},h::SVector{N,T},fit=ones(SVector{N,Bool}))

Description
"""
function refineCell!(dom::AbstractNestedDomain{N,T},uppos::SVector{N,Int},h::SVector{N,T},fit=ones(SVector{N,Bool})) where {N,T}
    isin = dom.table[uppos...]
    if isin == 0
        return
    elseif isin == -1
        return
    end
    limits = get_cell_from_pos(dom, uppos)
    dom.table[uppos...] = -1
    depth = dom.depth+1
    isroot = false
    isleaf = true
    updomain = Ref(dom)

    dims = Vector{Int}(undef,N) #dimensions of the table

    newh = Vector{T}(undef,N)
    newSize = Vector{T}(undef,N)
    for i = 1:N
        if fit[i] || periods[i][0] !=0
            dims[i] = round(Int, UT.size(limits)[i]/h[i])
            newh[i] = UT.size(limits)[i]/dims[i]
            newSize[i] = UT.size(limits)[i]
        else
            dims[i] = trunc(Int, UT.size(limits)[i]/h[i])
            newSize[i] = dims[i]*h[i]
            newh[i] = h[i]
        end
    end

    h = SVector{N}(newh)
    limits = UT.buildRectangle2(limits.lb, SVector{N}(newSize))

    table = ones(Int,dims...)
    subdomains = Dict{SVector{N,Int}, Ref}()
    # cells = Dict{SVector{N,Int}, Int}()

    dom.subdomains[uppos] = Ref(NestedDomain(limits,h,depth,isroot,isleaf,dom.isperiodic,dom.periods,updomain,uppos,table,subdomains))
    return enum_centers(dom.subdomains[uppos][]), [UT.center(limits)]
end

# CHECKED
# Assume that center is in the selected coordinates
"""
    addFreeDomain!(dom::AbstractNestedDomain{N,T}, uppos::SVector, center::SVector)

Description
"""
function addFreeDomain!(dom::AbstractNestedDomain{N,T}, uppos::SVector, center::SVector) where {N,T}
    added_centers = SVector{N,T}[]
    removed_centers = SVector{N,T}[]

    limits = get_cell_from_pos(dom, uppos)

    if all(center .≈ limits.lb) || all(center .≈ limits.ub)
        return added_centers, removed_centers
    end

    push!(removed_centers, UT.center(get_cell_from_pos(dom, uppos)))

    dom.table[uppos...] = -1
    depth = dom.depth+1
    isroot = false
    isleaf = true
    updomain = Ref(dom)

    dims = Vector{Int}(undef,N) #dimensions of the table

    for i = 1:N
        if center[i] ≈ limits.lb[i] || center[i] ≈ limits.ub[i]
            dims[i] = 1
        else
            dims[i] = 2
        end
    end

    table = ones(Int,dims...) #check syntax

    subdomains = Dict{SVector{N,Int}, Ref}() # empty since all cells are labeled as 1 or 0

    # cells = Dict{SVector{N,Int}, Int}() # empty for now. It can be constructed together with the automaton.
    #                                 # it could also be computed rn but i dunno if it's the best thing to do
    dom.subdomains[uppos] = Ref(FreeNestedDomain(limits,depth,isroot,isleaf,dom.isperiodic,dom.periods,updomain,uppos,center,table,subdomains))
    append!(added_centers, enum_centers(dom.subdomains[uppos][]))
    return added_centers, removed_centers
end
# WRONG
# TODO:Fix
"""
    addFreeDomain!(dom::AbstractNestedDomain{N,T}, uppos::SVector, R::UT.HyperRectangle)

Description
"""
function addFreeDomain!(dom::AbstractNestedDomain{N,T}, uppos::SVector, R::UT.HyperRectangle) where {N,T}
    added_centers = SVector{N,T}[]
    removed_centers = SVector{N,T}[]

    limits = get_cell_from_pos(dom, uppos)

    do_intersect, R_int = UT.computeIntersection!(limits, R)
    if uppos == SVector(2,3,1)
        println(limits)
        println(R_int)
    end
    if do_intersect == 0 || do_intersect == 2 || do_intersect == 4
        return added_centers, removed_centers
    end

    Mcenter = Vector{T}(undef,N)
    for i = 1 : N
        if !UT.leaq(R_int.lb[i], limits.lb[i])
            Mcenter[i] = R_int.lb[i]
        else
            Mcenter[i] = R_int.ub[i]
        end
    end
    center = SVector{N}(Mcenter)

    added, removed = addFreeDomain!(dom,uppos,center)
    append!(added_centers, added)
    append!(removed_centers, removed)

    if !isempty(added)#this means also that both is_lb_inside and is_ub_inside are true
        newsubdom = dom.subdomains[uppos][]
        newpos = SVector{N}(size(newsubdom.table))
        newcenter = R_int.ub
        added, removed = addFreeDomain!(newsubdom,newpos,newcenter)
        append!(added_centers, added)
        append!(removed_centers, removed)
    end

    return added_centers, removed_centers
end

# CHECKED
"""
    addObstacleDomain!(dom::AbstractNestedDomain{N,T}, uppos::SVector{N,Int}, obstacle::UT.HyperRectangle)

Description
"""
function addObstacleDomain!(dom::AbstractNestedDomain{N,T}, uppos::SVector{N,Int}, obstacle::UT.HyperRectangle) where {N,T}
    added_centers = SVector{N,T}[]
    removed_centers = SVector{N,T}[]

    limits = get_cell_from_pos(dom, uppos)
    do_intersect, intersection = UT.computeIntersection!(limits, obstacle)
    if do_intersect == 0
        return added_centers, removed_centers
    elseif do_intersect == 2 || do_intersect == 4 # simply remove cell
        dom.table[uppos...] = 0
        push!(removed_centers, UT.center(get_cell_from_pos(dom, uppos)))
        return added_centers, removed_centers
    end

    push!(removed_centers, UT.center(get_cell_from_pos(dom, uppos)))

    dom.table[uppos...] = -1
    depth = dom.depth+1
    isroot = false
    isleaf = true
    updomain = Ref(dom)

    dims = Vector{Int}(undef,N) #dimensions of the table
    obstacle_pos = Vector{Int}(undef,N)
    for i = 1:N
        is_on_lb = false
        is_on_ub = false
        if intersection.lb[i] ≈ limits.lb[i] #intersects the i-th lower border
            is_on_lb = true
        end
        if intersection.ub[i] ≈ limits.ub[i] #intersects the i-th upper border
            is_on_ub = true
        end
        if is_on_lb
            obstacle_pos[i] = 1 #is on the i-th lb
            if (is_on_ub) dims[i] = 1
            else dims[i] = 2 end
        else
            obstacle_pos[i] = 2 #is on the i-th ub or in the middle
            if (is_on_ub) dims[i] = 2
            else dims[i] = 3 end
        end
    end

    table = ones(Int,dims...)
    table[obstacle_pos...] = 0

    subdomains = Dict{SVector{N,Int}, Ref}() # empty since all cells are labeled as 1 or 0

    # cells = Dict{SVector{N,Int}, Int}() # empty for now. It can be constructed together with the automaton.
    #                                 # it could also be computed rn but i dunno if it's the best thing to do
    dom.subdomains[uppos] = Ref(ObstacleNestedDomain(limits,depth,isroot,isleaf,dom.isperiodic,dom.periods,updomain,uppos,intersection,SVector{N}(obstacle_pos),table,subdomains))
    append!(added_centers, enum_centers(dom.subdomains[uppos][]))
    return added_centers, removed_centers
end

# CHECKED
"""
    removeObstacle!(dom::AbstractNestedDomain{N,T},obstacle::UT.HyperRectangle,bringinperiod=true)

Description
"""
function removeObstacle!(dom::AbstractNestedDomain{N,T},obstacle::UT.HyperRectangle,bringinperiod=true) where {N,T}
    if bringinperiod
        obstacles = bring_into_period(dom,obstacle)
        return removeObstacles!(dom,obstacles,false)
    end

    added_centers = SVector{N,T}[]
    removed_centers = SVector{N,T}[]

    does_int, intersection = UT.computeIntersection!(dom.limits, obstacle)
    if does_int == 0
        return added_centers, removed_centers
    elseif does_int == 2 || does_int == 4
        removeCoordFromDomain!(dom.updomain[], dom.updomain_pos) #remove domain from upper one
                        # i haven't removed the reference to the subdomain, updated the mdp etc.
        push!(removed_centers, UT.center(get_cell_from_pos(dom.updomain[], dom.updomain_pos)))
        return added_centers, removed_centers
    end

    # corners = UT.get_corners(intersection)
    _, obst_lb_pos = get_pos_from_coord(dom,intersection.lb,false,false)
    _, obst_ub_pos = get_pos_from_coord(dom,intersection.ub,false,false, get_lower_one=true)

    # remove all the cells
    t_size = size(dom.table)
    obstacle_size = obst_ub_pos .- obst_lb_pos .+1
    M = prod(obstacle_size)

    curr_Mpos = Vector{Int}(undef,N) #current position of the cell that i want to "remove"
    for k = 0 : M-1 # this loop should allow to visit all the cells in the intersection
        tmp = k
        for i = 1 : N # maps k into a pos
            curr_Mpos[i] = tmp % obstacle_size[i] + obst_lb_pos[i]
            if curr_Mpos[i] > t_size[i]
                curr_Mpos[i] = t_size[i]
            end
            tmp = trunc(Int, tmp / obstacle_size[i])
        end
        curr_pos = SVector{N}(curr_Mpos)
        isin = dom.table[curr_pos...]
        if isin == -1
            added, removed = removeObstacle!(dom.subdomains[curr_pos][],intersection,false)
            append!(added_centers, added)
            append!(removed_centers, removed)
        elseif isin == 0 #do nothing
        elseif !any(UT.leaq.(curr_pos,obst_lb_pos)) && !any(UT.greaq.(curr_pos, obst_ub_pos)) #the cell is completely inside the obstacle
            removeCoordFromDomain!(dom,curr_pos)
            push!(removed_centers, UT.center(get_cell_from_pos(dom, curr_pos)))
        else
            added, removed = addObstacleDomain!(dom, curr_pos, intersection)
            append!(added_centers, added)
            append!(removed_centers, removed)
        end
    end
    return added_centers, removed_centers
end
# # CHECKED
# function removeObstacle!(dom::ObstacleNestedDomain{N,T},obstacle::UT.HyperRectangle,bringinperiod=true) where {N,T}
#     if bringinperiod
#         obstacles = bring_into_period(dom,obstacle)
#         removeObstacles!(dom,obstacles,false)
#         return
#     end
#
#     does_int, intersection = UT.computeIntersection!(dom.limits, obstacle)
#     if does_int == 0
#         return
#     elseif does_int == 2 || does_int == 4
#         removeCoordFromDomain!(dom.updomain[], dom.updomain_pos) #remove domain from upper one
#                         # i haven't removed the reference to the subdomain, updated the mdp etc.
#         return
#     end
#
#     _, obst_lb_pos = get_pos_from_coord(dom,intersection.lb,false,false)
#     _, obst_ub_pos = get_pos_from_coord(dom,intersection.ub,false,false)
#
#     t_size = size(dom.table)
#     M = prod(t_size)
#
#     curr_Mpos = Vector{Int}(undef,N) #current position of the cell that i want to operate on
#     for k = 0 : M-1 # this loop should allow to visit all the cells in the intersection
#         tmp = k
#         for i = 1 : N # maps k into a pos
#             curr_Mpos[i] = tmp % t_size[i] + 1
#             tmp = trunc(Int, tmp / t_size[i]) + 1
#         end
#         curr_pos = SVector{N}(curr_Mpos)
#         isin = dom.table[curr_pos...]
#         if isin == -1
#             removeObstacle!(dom.subdomains[curr_pos][],intersection)
#         elseif isin == 0 #do nothing
#         elseif all(curr_pos .> obst_lb_pos) && all(curr_pos .< obst_ub_pos) #the cell is completely inside the obstacle
#             removeCoordFromDomain!(dom,curr_pos)
#         else
#             addObstacleDomain!(dom, curr_pos, intersection)
#         end
#     end
# end
# CHECKED
"""
    removeObstacles!(dom::AbstractNestedDomain{N,T},obstacles,bringinperiod=true)

Description
"""
function removeObstacles!(dom::AbstractNestedDomain{N,T},obstacles,bringinperiod=true) where {N,T}
    added_centers = SVector{N,T}[]
    removed_centers = SVector{N,T}[]
    for obstacle in obstacles
        added,removed = removeObstacle!(dom,obstacle,bringinperiod)
        append!(added_centers,added)
        append!(removed_centers,removed)
    end
    return added_centers,removed_centers
end

# CHECKED
# return the length of the period along dimension i
"""
    period_length(periods::SVector, i)

Description
"""
function period_length(periods::SVector, i)
    return periods[i][2]
end
"""
    period_length(dom::AbstractNestedDomain, i)

Description
"""
function period_length(dom::AbstractNestedDomain, i)
    return dom.periods[i][2]
end
"""
    period_lb(periods::SVector, i)

Description
"""
function period_lb(periods::SVector, i)
    return periods[i][1]
end
"""
    period_lb(dom::AbstractNestedDomain, i)

Description
"""
function period_lb(dom::AbstractNestedDomain, i)
    return dom.periods[i][1]
end
"""
    period_ub(periods::SVector, i)

Description
"""
function period_ub(periods::SVector, i)
    return periods[i][1]+periods[i][2]
end
"""
    period_ub(dom::AbstractNestedDomain, i)

Description
"""
function period_ub(dom::AbstractNestedDomain, i)
    return dom.periods[i][1]+dom.periods[i][2]
end
"""
    is_periodic(periods::SVector, i)

Description
"""
function is_periodic(periods::SVector, i)
    return periods[i][2] > 0
end
"""
    is_periodic(dom::AbstractNestedDomain, i)

Description
"""
function is_periodic(dom::AbstractNestedDomain, i)
    return dom.periods[i][2] > 0
end

# CHECKED
# helper function (but it can be useful)
function _shift_into_period(period::SVector{2,T}, xi::T) where T
    if period[2] == 0 || period[1]<xi<period[1]+period[2]
        return xi
    end
    return period[1] + ( (xi-period[1]) % period[2] )
end
# CHECKED
# helper function (but it can be useful)
"""
    bring_into_period(periods, x::SVector{N,T})

Description
"""
function bring_into_period(periods, x::SVector{N,T}) where {N,T}
    m = Vector{T}(undef,N)
    for i = 1 : N
        m[i] = _shift_into_period(periods[i], x[i])
    end
    return SVector{N}(m)
end
# CHECKED
# returns the point shifted in the period
"""
    bring_into_period(dom::AbstractNestedDomain, x::SVector{N,T})

Description
"""
function bring_into_period(dom::AbstractNestedDomain, x::SVector{N,T}) where {N,T}
    if !dom.isperiodic
        return x
    end
    return bring_into_period(dom.periods,x)
end

# CHECKED
function _split_one_dimension(lb_i,ub_i,period)
    #helper function
    if UT.greaq(ub_i-lb_i, period[2])
        return [SVector(period[1],period[1]+period[2])]
    else
        lb_i = _shift_into_period(period, lb_i)
        ub_i = _shift_into_period(period, ub_i)
        if UT.leaq(lb_i,ub_i)
            return [SVector(lb_i,ub_i)]
        else
            return [SVector(period[1],ub_i),SVector(lb_i,period[1]+period[2])]
        end
    end
end
# CHECKED
function _recursive_split!(L,periods::SVector{N,SVector{2,T}},rect::UT.HyperRectangle,curr_lb::Vector{T},curr_ub::Vector{T},i::Int) where {N,T}
    #helper function
    if i > N #i is the current dimension that i'm splitting
        push!(L, UT.HyperRectangle(SVector{N}(curr_lb), SVector{N}(curr_ub)))
        return
    end
    new_lb = copy(curr_lb)
    new_ub = copy(curr_ub)
    if !is_periodic(periods,i)
        new_lb[i] = rect.lb[i]
        new_ub[i] = rect.ub[i]
        _recursive_split!(L,periods,rect,new_lb,new_ub,i+1)
    else
        intervals = _split_one_dimension(rect.lb[i],rect.ub[i],periods[i])
        for interval in intervals
            new_lb[i] = interval[1]
            new_ub[i] = interval[2]
            _recursive_split!(L,periods,rect,new_lb,new_ub,i+1)
        end
    end
end
# CHECKED
# returns a list of rects which are the result of shifting the input rect in the period
# TODO: it can be improved so that ""it does not to split each dimension 2^N times"" and
#   +make it non-recursive
#   +precompute the dimension of L
"""
    bring_into_period(periods::SVector{N,SVector{2,T}}, rect::UT.HyperRectangle)

Description
"""
function bring_into_period(periods::SVector{N,SVector{2,T}}, rect::UT.HyperRectangle) where {N,T}
    L = typeof(rect)[]
    curr_lb = Vector{T}(undef,N)
    curr_ub = Vector{T}(undef,N)
    _recursive_split!(L,periods,rect,curr_lb,curr_ub,1)
    return L
end
"""
    bring_into_period(dom::AbstractNestedDomain{N,T}, rect::UT.HyperRectangle)

Description
"""
function bring_into_period(dom::AbstractNestedDomain{N,T}, rect::UT.HyperRectangle) where {N,T}
    if dom.isperiodic
        return bring_into_period(dom.periods, rect)
    else
        return [rect]
    end
end

# CHECKED
function _shape_rect2D(rect::UT.HyperRectangle)
    Plots.Shape([(rect.lb[1],rect.lb[2]),(rect.lb[1],rect.ub[2]),(rect.ub[1],rect.ub[2]),(rect.ub[1],rect.lb[2])])
end
# CHECKED
"""
    plot_domain!(dom::AbstractNestedDomain{N,T},fig,until_depth)

Description
"""
function plot_domain!(dom::AbstractNestedDomain{N,T},fig,until_depth) where {N,T}#2D
    dims = size(dom.table)
    pos = ones(Int,N)
    for i = 1 : dims[1]
        pos[1] = i
        for j = 1 : dims[2]
            pos[2] = j
            iscell = dom.table[pos...]
            if iscell == 1
                cell = get_cell_from_pos(dom, SVector{N}(pos))
                Plots.plot!(_shape_rect2D(cell), opacity=0.2,color=:green)
            elseif iscell == 0
                cell = get_cell_from_pos(dom, SVector{N}(pos))
                Plots.plot!(_shape_rect2D(cell), opacity=0.2,color=:red)
            else
                if until_depth == dom.depth
                    cell = get_cell_from_pos(dom, SVector{N}(pos))
                    Plots.plot!(_shape_rect2D(cell), opacity=0.2,color=:yellow)
                else
                    subdom = dom.subdomains[pos][]
                    plot_domain!(subdom,fig,until_depth)
                end
            end
        end
    end
    return
end
# CHECKED
"""
    plot_domain(dom::AbstractNestedDomain, until_depth::Int=-1; init=nothing, target=nothing, fig=Plots.plot(aspect_ratio = 1,legend = false))

Description
"""
function plot_domain(dom::AbstractNestedDomain, until_depth::Int=-1; init=nothing, target=nothing, fig=Plots.plot(aspect_ratio = 1,legend = false))
    plot_domain!(dom,fig,until_depth)
    if init != nothing
        Plots.plot!(_shape_rect2D(init), opacity=1,color=:black)
    end
    if target != nothing
        Plots.plot!(_shape_rect2D(target), opacity=0.8,color=:blue)
    end
    display(fig)
    return
end

"""
    addSet!(dom::AbstractNestedDomain{N,T},R::UT.HyperRectangle,bringinperiod=true)

Description
"""
function addSet!(dom::AbstractNestedDomain{N,T},R::UT.HyperRectangle,bringinperiod=true) where {N,T}
    if bringinperiod
        Rlist = bring_into_period(dom,R)
        return addSets!(dom,Rlist,false)
    end

    added_centers = SVector{N,T}[]
    removed_centers = SVector{N,T}[]

    does_int, R_int = UT.computeIntersection!(dom.limits, R)
    if does_int == 0 || does_int == 2 || does_int == 4 # limits <: R or doesnt int
        return added_centers, removed_centers
    end

    _, R_lb_pos = get_pos_from_coord(dom,R_int.lb,false,false)
    _, R_ub_pos = get_pos_from_coord(dom,R_int.ub,false,false, get_lower_one=true)

    t_size = size(dom.table)
    R_size = R_ub_pos .- R_lb_pos .+1
    M = prod(R_size)

    curr_Mpos = Vector{Int}(undef,N)
    for k = 0 : M-1 # this loop should allow to visit all the cells in the intersection
        tmp = k
        for i = 1 : N # maps k into a pos
            curr_Mpos[i] = tmp % R_size[i] + R_lb_pos[i]
            if UT.greaq(curr_Mpos[i], t_size[i])
                curr_Mpos[i] = t_size[i]
            end
            tmp = trunc(Int, tmp / R_size[i])
        end
        curr_pos = SVector{N}(curr_Mpos)
        isin = dom.table[curr_pos...]
        if isin == -1
            added, removed = addSet!(dom.subdomains[curr_pos][],R_int,false)
            append!(added_centers, added)
            append!(removed_centers, removed)
        elseif isin == 0 #do nothing
        elseif !any(UT.leaq.(curr_pos,R_lb_pos)) && !any(UT.greaq.(curr_pos, R_ub_pos)) #the cell is completely inside the set
            #do nothing
        else
            added, removed = addFreeDomain!(dom, curr_pos, R_int)
            append!(added_centers, added)
            append!(removed_centers, removed)
        end
    end
    return added_centers, removed_centers
end
"""
    addSets!(dom::AbstractNestedDomain{N,T},Rlist,bringinperiod=true)

Description
"""
function addSets!(dom::AbstractNestedDomain{N,T},Rlist,bringinperiod=true) where {N,T}
    added_centers = SVector{N,T}[]
    removed_centers = SVector{N,T}[]
    for R in Rlist
        added,removed = addSet!(dom,R,bringinperiod)
        append!(added_centers,added)
        append!(removed_centers,removed)
    end
    return added_centers,removed_centers
end

function _rec_enum_centers_in_rect!(dom::AbstractNestedDomain{N,T},X,r) where {N,T}

    _, r_lb_pos = get_pos_from_coord(dom,r.lb,false,false)
    _, r_ub_pos = get_pos_from_coord(dom,r.ub,false,false, get_lower_one=true)

    t_size = size(dom.table)
    r_size = r_ub_pos .- r_lb_pos .+1
    M = prod(r_size)

    curr_Mpos = Vector{Int}(undef,N) #current position of the cell
    for k = 0 : M-1 # this loop should allow to visit all the cells in the rect
        tmp = k
        for i = 1 : N # maps k into a pos
            curr_Mpos[i] = tmp % r_size[i] + r_lb_pos[i]
            if UT.greaq(curr_Mpos[i], t_size[i])
                curr_Mpos[i] = t_size[i]
            end
            tmp = trunc(Int, tmp / r_size[i])
        end

        curr_pos = SVector{N}(curr_Mpos)
        isin = dom.table[curr_pos...]
        if isin == 1
            cell = get_cell_from_pos(dom,curr_pos)
            push!(X,UT.center(cell))
        elseif isin == -1
            cell = get_cell_from_pos(dom,curr_pos)
            _rec_enum_centers_in_rect!(dom.subdomains[curr_pos][],X,cell)
        end
    end
    return
end
# NOTE: does not take into account of periodicity
"""
    enum_centers_in_rect(dom::AbstractNestedDomain{N,T},r::UT.HyperRectangle)

Description
"""
function enum_centers_in_rect(dom::AbstractNestedDomain{N,T},r::UT.HyperRectangle) where {N,T}
    does_int, newr = UT.computeIntersection!(dom.limits, r)

    X = SVector{N,T}[]
    if does_int != 0
        _rec_enum_centers_in_rect!(dom,X,newr)
    end
    return X
end
"""
    enum_centers(dom::AbstractNestedDomain)

Description
"""
function enum_centers(dom::AbstractNestedDomain)
    return enum_centers_in_rect(dom, dom.limits)
end

#assume that they do not overlap
"""
    enum_centers_in_rects(dom::AbstractNestedDomain{N,T},R_list)

Description
"""
function enum_centers_in_rects(dom::AbstractNestedDomain{N,T},R_list) where {N,T}
    X = SVector{N,T}[]
    for r in R_list
        append!(X,enum_centers_in_rect(dom,r))
    end
    return X
end
