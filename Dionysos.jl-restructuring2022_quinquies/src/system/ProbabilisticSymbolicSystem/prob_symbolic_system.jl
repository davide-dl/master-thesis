abstract type AbstractProbSymbolicSys{Nx,Nu,T<:Int} <: AbstractSystem{Nx,Nu,T} end
abstract type AbstractNode end
abstract type AbstractAction end
abstract type AbstractTransition end

struct Transition <: AbstractTransition
    prob::Float64
    # other stats that can be useful?
end
struct Action <: AbstractAction
    # symbol::Int
    # is_allowed::Bool
    trans::Dict{Int,Transition} # key = symbol of the final_state
    # other stats that can be useful eg. q-value
end
# action with time_step
struct ActionT <: AbstractAction
    trans::Dict{Int,Transition} # key = symbol of the final_state
    tstep::Real
end
struct Node <: AbstractNode
    # symbol::Int
    # cell::UT.HyperRectangle #or something equivalent (eg. pointer to cell)
    # n_actions::Int
    actions::Dict{Int, AbstractAction} # key = symbol of the input
    entering_trans::Dict{Int,Set{Int}} # entering transitions
                               # key = initial state
                               # value = actions
    # flag::Int16
    # other stats that can be useful eg. score
end
struct DictProbAutomaton{Nx,Nu,T} <: AbstractProbSymbolicSys{Nx,Nu,T}
    nodes::Dict{Int,Node} # key = symbol of the initial_state
    # n_nodes::Int #number of nodes
    #t_step ?
end

struct DictProbAutomatonT{Nx,Nu,T} <: AbstractProbSymbolicSys{Nx,Nu,T}
    nodes::Dict{Int,Node} # key = symbol of the initial_state
    # n_nodes::Int #number of nodes
    #t_step ?
end

function buildEmptyDictProbAutomaton()
    nodes = Dict{Int,Node}()
    return DictProbAutomaton{1,1,Int}(nodes)
end
function buildEmptyNode()
    actions = Dict{Int,AbstractAction}()
    entering_trans = Dict{Int,Vector{Int}}()
    return Node(actions,entering_trans)
end
function buildEmptyAction()
    transitions = Dict{Int,Transition}()
    return Action(transitions)
end
function buildEmptyActionT(tstep)
    transitions = Dict{Int,Transition}()
    return ActionT(transitions,tstep)
end
function buildTransition(prob)
    return Transition(prob)
end

function addEnteringTrans!(sys::DictProbAutomaton, s1::Int, a::Int, s2::Int)
    # if get(sys.nodes,s2,false) == false
    #     sys.nodes[s2] = buildEmptyNode()
    # end
    curr_trans = get(sys.nodes[s2].entering_trans, s1, Set{Int}())
    sys.nodes[s2].entering_trans[s1] = push!(curr_trans, a)
end

function get_prob(sys::DictProbAutomaton,s1::Int,a::Int,s2::Int)
    return sys.nodes[s1].actions[a].trans[s2].prob
end
function list_trans(sys::DictProbAutomaton,s1::Int,a::Int)
    return [s2 for (s2,_) in sys.nodes[s1].actions[a].trans]
end
function list_actions(sys::DictProbAutomaton,s1::Int)
    return [a for (a,_) in sys.nodes[s1].actions]
end
function list_nodes(sys::DictProbAutomaton)
    return [s for (s,_) in sys.nodes]
end
function list_entering_nodes(sys::DictProbAutomaton,s2::Int)
    return [s1 for (s1,_) in sys.nodes[s2].entering_trans]
end
