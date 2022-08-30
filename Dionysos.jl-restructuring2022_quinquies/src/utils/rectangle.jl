struct HyperRectangle{VT}
    lb::VT
    ub::VT
end

# Changed all(f(x[i]) for i in eachindex(x)) to all(f.(x))
# See test_performances.
function Base.in(x, rect::HyperRectangle)
    return all(leaq.(rect.lb, x)) && all(leaq.(x, rect.ub))
end

function Base.isempty(rect::HyperRectangle)
    return !all(leaq.(rect.lb, rect.ub))
end

function Base.intersect(a::HyperRectangle, b::HyperRectangle)
    return HyperRectangle(max.(a.lb, b.lb), min.(a.ub, b.ub))
end

function Base.issubset(a::HyperRectangle, b::HyperRectangle)
    return all(greaq.(a.lb, b.lb)) && all(leaq.(a.ub, b.ub))
end

function Base.isequal(a::HyperRectangle, b::HyperRectangle)
    return (all(a.lb .≈ b.lb) && all(a.ub .≈ b.ub))
end


function size(rect::HyperRectangle)
    return rect.ub.-rect.lb
end
function center(rect::HyperRectangle)
    return rect.lb .+ size(rect)./2
end
function dim(rect::HyperRectangle)
    return length(rect.lb)
end
function buildRectangle(lb,ub)
    return HyperRectangle(lb,ub)
end
function buildRectangle2(lb,size)
    return HyperRectangle(lb,lb.+size)
end

function isok(r::HyperRectangle)
    return all(greaq.(r.ub, r.lb))
end

# returns 0 if they do not intersect or if the instersection has volume=0
#         1 if they intersect
#         2 if they are equal
#         3 if b is contained in a
#         4 if a is contained in b
leaq(a,b) = (a <= b) || (a ≈ b)
greaq(a,b) = (a >= b) || (a ≈ b)
# leaq(a,b) = (a <= b) || isapprox(a,b,atol=1e-8)
# greaq(a,b) = (a >= b) || isapprox(a,b,atol=1e-8)
function computeIntersection!(a::HyperRectangle,b::HyperRectangle)
    output = intersect(a, b)
    if any(greaq.(output.lb, output.ub)) #isempty
        return 0, output
    end
    if all(a.lb ≈ b.lb) && all(a.ub ≈ b.ub)
        return 2, output
    end
    lba = false
    lbb = false
    uba = false
    ubb = false
    if all(leaq.(a.lb, b.lb))
        lba = true
    end
    if all(greaq.(a.lb, b.lb))
        lbb = true
    end
    if all(leaq.(a.ub, b.ub))
        uba = true
    end
    if all(greaq.(a.ub, b.ub))
        ubb = true
    end
    if lba & ubb
        return 3, output
    elseif lbb & uba
        return 4, output
    else
        return 1, output
    end
end

# returns vector containing all the corners/vertices of the hyperrectangle
function get_corners(rect::HyperRectangle{SVector{N,T}}) where {N,T}
    pointVec = Vector{SVector{N,T}}(undef,2^N)
    p = MVector{N,T}(undef)
    lwb = rect.lb
    upb = rect.ub
    for i = 1 : 2^N
        bits  = bitstring(Int16(i-1))
        for j = 1 : N
            p[j] = (bits[16+1-j] == '0') ? lwb[j] : upb[j]
        end
        pointVec[i] = SVector(p)
    end
    return pointVec
end
