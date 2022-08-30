struct ListDomain{N,T} <: DomainType{N,T}
    inputs::Vector
end

function enum_coords(d::ListDomain)
    return d.inputs
end

function build1DListDomain(inputs::Vector{SVector{N,T}}) where {N,T}
    return ListDomain{N,T}(inputs)
end
# function build1DListDomain(inputs)
#     T = typeof(inputs[1][1])
#     return ListDomain{1,T}(inputs)
# end
