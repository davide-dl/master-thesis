#NOTE: I dunno if this is a fast implementation

abstract type AbstractSymbolicDomain{N,T} <: DomainType{N,T} end # It could be a "normal" nested domain or a FreeNestedDomain

# each index of the vector is 0 if the symbol is NOT in the domain
# otherwise to each element can be assigned a Int value used as flag
# 1 is the default flag
# -1 = prohibited
# -2 = unreachable from init
# free_indexes stores the elems of the
struct SimpleSymbolicDomain{N,T<:Int} <: AbstractSymbolicDomain{N,T}
    elems::Vector{Int16}
    free_indexes::BitSet
    # last_index::Int64
end

function buildSymbolicDomain()
    return SimpleSymbolicDomain{1,Int64}(Vector{Bool}(),BitSet())
end

function buildSymbolicDomain(num_elems::Int)
    return SimpleSymbolicDomain{1,Int64}(Vector{Bool}(trues(num_elems)),BitSet())
end

function add_next_s!(d::SimpleSymbolicDomain)
    if isempty(d.free_indexes)
        push!(d.elems,1)
        return lastindex(d.elems)
    end
    return pop!(d.free_indexes)
end
function remove_s!(d::SimpleSymbolicDomain, s::Int)
    if get(d.elems,s,0) != 0
        d.elems[s] = 0
        push!(d.free_indexes,s)
    end
end
function get_next_s(d::SimpleSymbolicDomain)
    if isempty(d.free_indexes)
        return lastindex(d.elems)+1
    end
    return first(d.free_indexes)
end
function get_last_s(d::SimpleSymbolicDomain)
    return findlast(d.elems)
end

function set_s!(d::SimpleSymbolicDomain, s::Int, flag::Int)
    if flag == 0
        return remove_s!(d,s)
    end
    if get(d.elems,s,0) == 0
        delete!(d.free_indexes,s)
    end
    d.elems[s] = flag
end


function enum_elems(d::SimpleSymbolicDomain)
    return [i for (i,e) in enumerate(d.elems) if e!=0]
end
function enum_elems(d::SimpleSymbolicDomain, flag)
    return [i for (i,e) in enumerate(d.elems) if e==flag]
end
function enum_elems(d::SimpleSymbolicDomain, condition::Function)
    return [i for (i,e) in enumerate(d.elems) if condition(e)]
end

function is_in(d::SimpleSymbolicDomain, s::Int)
    if d.elems[s] != 0
        return true
    else
        return false
    end
end

# function add!(d::SimpleSymbolicDomain,s::Int)
#     d.elems[s] = true
#     if lengthfree_indexes
# end
#
# function add_next!(d::SimpleSymbolicDomain)
#     push!(d.elems,s)
# end
#
# function in(s::Int, d::SimpleSymbolicDomain)
#     return in(s,d.elems)
# end
