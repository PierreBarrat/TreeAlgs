@kwdef struct FitchState{T}
    state::Set{T}
end

FitchState(::Val{T}) where T = FitchState(Set{T}())
FitchState(i::Int, a::T) where T = FitchState(Set(a))
Base.isempty(x::FitchState) = isempty(x.state)
Base.length(x::FitchState) = length(x.state)

@kwdef mutable struct FitchData{T} <: TreeTools.TreeNodeData
    observed_sequence::Vector{T}
    reconstructed_sequence::Vector{T}
    state::FitchState{T}
end

FitchData(::Val{T}) where T = FitchData(T[], T[], FitchState(Val{T}()))
FitchData(os::BioSequence, rs::BioSequence) = FitchData(collect(os), collect(rs))
FitchData(os::T, rs::T) where T = FitchData(os, rs, FitchState(Val{eltype(T)}()))
function FitchData(observed_sequence :: AbstractString)
    return FitchData(collect(observed_sequence), Vector{Char}(undef, 0))
end
function FitchData(observed_sequence :: AbstractVector{T}) where T
    return FitchData(observed_sequence, Vector{T}(undef, 0))
end
FitchData(os::BioSequence) = FitchData(collect(os))

function Base.isempty(x::FitchData)
    return isempty(x.observed_sequence) &&
        isempty(x.reconstructed_sequence) &&
        isempty(x.state)
end

let pos::Int = 0
    global current_pos() = pos
    global reset_pos!() = (pos = 0)
    global set_pos!(i) = (pos = i)
end

"""
    init_fitchstate!(dat::FitchData)

Initialize the fitch state at `current_pos()` using the observed sequence.
If `dat.observed_sequence` is empty, the state remains empty.
If it is not empty but shorter than `current_pos()`, an error is thrown.
"""
function init_fitchstate!(dat::FitchData)
    if !isempty(dat.observed_sequence)
        if length(dat.observed_sequence) > current_pos()
            throw(ErrorException("
                Sequence shorter than expected:\
                wanted at least $(current_pos()), but $(length(dat.observed_sequence))
            "))
        end
        dat.state = Set(dat.observed_sequence[current_pos()])
    end
    return dat.state
end
init_fitchstate!(n::TreeNode{<:FitchData}) = init_fitchstate!(data(n))

function fitch_up!(n::TreeNode{<:FitchData})
    # recursive call
    for c in children(n)
        fitch_up!(c)
    end
    # actual computation: build the state of n using states of children
    pull_downstream_state!(n)
end
function pull_downstream_state!(n::TreeNode)
    # trying intersect
    # union with first child (to init.) then intersect with the rest
    c1, other_children = firstrest(children(n))
    union!(data(n).state.state, c1.data.state.state)
    for c in other_children
        intersect!(data(n).state.state, data(c).state.state)
    end
    # if intersect is empty, go for union!
    if isempty(data(n).state.state)
        for c in children(n)
            union!(data(n).state.state, data(c).state.state)
        end
    end

    return nothing
end

function sample_from_sate!(dat::FitchData)

end



