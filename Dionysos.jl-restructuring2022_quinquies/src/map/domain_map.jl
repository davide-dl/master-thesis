"""
    DictDomainMap{N1,T1,N2,T2} <: AbstractDomainMap{N1,T1,N2,T2}

Description
"""
struct DictDomainMap{N1,T1,N2,T2} <: AbstractDomainMap{N1,T1,N2,T2}
    dom1::DO.DomainType{N1,T1}
    dom2::DO.DomainType{N2,T2}
    dict1::Dict #1->2
    dict2::Dict #2->1
end
struct StaticDictMap{N1,T1,N2,T2} <: AbstractDomainMap{N1,T1,N2,T2}
    dom1::DO.DomainType{N1,T1}
    dom2::DO.DomainType{N2,T2}
    dict1::Dict #1->2
    dict2::Dict #2->1
end
struct Nested2SymbolicMap{N1,T1,N2,T2} <: AbstractDomainMap{N1,T1,N2,T2}
    dom1::DO.DomainType{N1,T1}
    dom2::DO.DomainType{N2,T2}
    dict1::Dict{SVector,Int} #1->2
    dict2::Dict{Int,SVector}# 2->1: keys are symbols and values are (references to the domain, coordition)
end
struct Nested2SymbolicMap_2{N1,T1,N2,T2} <: AbstractDomainMap{N1,T1,N2,T2}
    dom1::DO.DomainType{N1,T1}
    dom2::DO.DomainType{N2,T2}
    dict1::Dict{SVector,Int} #1->2
    dict2::Dict{Int, Tuple{Ref,SVector}} # 2->1: keys are symbols and values are (references to the domain, coordition)
end

function map1(m::StaticDictMap,e1)
    return m.dict1[e1]
end
function map2(m::StaticDictMap,e2)
    return m.dict2[e2]
end

function map1(m::Nested2SymbolicMap,e1)
    return get(m.dict1,DO.get_cell_center_from_coord(m.dom1,e1), 0)
end
function map2(m::Nested2SymbolicMap,e2)
    return get(m.dict2,e2, nothing)
end

function map1(m::Nested2SymbolicMap{N1,T1,N2,T2},r::UT.HyperRectangle) where {N1,T1,N2,T2}
    centers = DO.enum_centers_in_rect(m.dom1,r)
    set = Int[]
    for c in centers
        push!(set, get(m.dict1,c, 0))
    end
    return set
end

function buildStaticMap(dom1::DO.DomainType{N1,T1}, dom2::DO.DomainType{N2,T2}) where {N1,T1,N2,T2} #I assume n_cells(dom1) == n_elems(dom2)
    elems1 = DO.enum_coords(dom1)
    elems2 = DO.enum_elems(dom2)
    dict1 = Dict(elems1 .=> elems2)
    dict2 = Dict(elems2 .=> elems1)
    return StaticDictMap{N1,T1,N2,T2}(dom1,dom2,dict1,dict2)
end
function buildStaticSymbolicMap(dom1::DO.DomainType{N1,T1}) where {N1,T1} #I assume n_cells(dom1) == n_elems(dom2)
    elems1 = DO.enum_coords(dom1)
    dom2 = DO.buildSymbolicDomain(length(elems1))
    elems2 = DO.enum_elems(dom2)
    dict1 = Dict(elems1 .=> elems2)
    dict2 = Dict(elems2 .=> elems1)
    return (StaticDictMap{N1,T1,1,Int64}(dom1,dom2,dict1,dict2), dom2)
end

function buildEmptyNested2SymbolicMap_2(dom1::DO.AbstractNestedDomain{N1,T1}, dom2::DO.AbstractSymbolicDomain{N2,T2}) where {N1,T1,N2,T2}
    dict1 = Dict{SVector,Int}() # center of the cell 2 symbol
    dict2 = Dict{Int, Tuple{Ref,SVector}}() # symbol to cell = (subdomain::Ref,coord::SVector{N,Int})
    return Nested2SymbolicMap_2{N1,T1,N2,T2}(dom1,dom2,dict1,dict2)
end
function buildNested2SymbolicMap_2(dom1::DO.AbstractNestedDomain{N1,T1}) where {N1,T1}
    X,C = DO.enum_centers_and_cells(dom1)
    dom2 = DO.buildSymbolicDomain(length(X))
    elems2 = DO.enum_elems(dom2)
    dict1 = Dict(X .=> elems2)
    dict2 = Dict(elems2 .=> C)
    return (Nested2SymbolicMap_2{N1,T1,1,Int64}(dom1,dom2,dict1,dict2), dom2)
end

function buildEmptyNested2SymbolicMap(dom1::DO.AbstractNestedDomain{N1,T1}, dom2::DO.AbstractSymbolicDomain{N2,T2}) where {N1,T1,N2,T2}
    dict1 = Dict{SVector,Int}() # center of the cell 2 symbol
    dict2 = Dict{Int,SVector}() # symbol to cell = (subdomain::Ref,coord::SVector{N,Int})
    return Nested2SymbolicMap{N1,T1,N2,T2}(dom1,dom2,dict1,dict2)
end
function buildNested2SymbolicMap(dom1::DO.AbstractNestedDomain{N1,T1}) where {N1,T1}
    X = DO.enum_centers(dom1)
    dom2 = DO.buildSymbolicDomain(length(X))
    elems2 = DO.enum_elems(dom2)
    dict1 = Dict(X .=> elems2)
    dict2 = Dict(elems2 .=> X)
    return (Nested2SymbolicMap{N1,T1,1,Int64}(dom1,dom2,dict1,dict2), dom2)
end
