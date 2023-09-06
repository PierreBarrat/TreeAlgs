mutable struct FitchData{T} <: TreeTools.TreeNodeData
    observed_sequence::Vector{T}
    reconstructed_sequence::Vector{T}
    state::Set{T}
end

function FitchData{T}() where T
    return FitchData{T}(Vector{T}(undef, 0), Vector{T}(undef, 0), Set{T}())
end

FitchData(::Val{T}) where T = FitchData(T[], T[], Set{T}())
FitchData(os::BioSequence, rs::BioSequence) = FitchData(collect(os), collect(rs))
FitchData(os::T, rs::T) where T = FitchData(os, rs, Set{eltype(T)}())
function FitchData(observed_sequence::AbstractString)
    return FitchData(collect(observed_sequence), Vector{Char}(undef, 0))
end
function FitchData(observed_sequence::AbstractVector{T}) where T
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

# """
#     fitch(tree::Tree, sequences::AbstractDict)
#     fitch(tree::tree, fastafile::AbstractString)
# """
function fitch(tree::Tree, sequences::AbstractDict)
    tree_fitch = sequences_to_tree(tree, sequences)
    fitch!(tree_fitch)
    return tree_fitch
end
function fitch(tree::Tree, fastafile::AbstractString)
    tree_fitch = fasta_to_tree(tree, fastafile)
    fitch!(tree_fitch)
    return tree_fitch
end
function fitch(tree::Tree{<:FitchData})
    tc = copy(tree)
    reset_fitch_states!(tc)
    fitch!(tc)

    return tc
end

function fitch!(tree::Tree{<:FitchData})
    check_tree(tree)
    L = let
        Ls = map(n -> length(n.data.observed_sequence), leaves(tree)) |> unique
        length(Ls) > 1 && error("All leaf sequences do not have the same length")
        first(Ls)
    end

    for pos in 1:L
        set_pos!(pos)
        foreach(init_fitchstate!, nodes(tree))
        fitch_up!(tree.root)
        fitch_down!(tree.root)
    end

    return nothing
end

function reset_fitch_states!(t::Tree{<:FitchData})
    return foreach(nodes(t)) do n
        n.data = FitchData(n.data.observed_sequence)
    end
end

#=
Check that tree is in a correct state to run Fitch
=#
function check_tree(tree)
    # all leaves should have an observed sequence
    for n in leaves(tree)
        if isempty(n.data.observed_sequence)
            throw(ErrorException("Leaf $(n.label) does not have an observed sequence"))
        end
    end
    # internals should not have an observed sequence
    for n in internals(tree)
        if !isempty(n.data.observed_sequence)
            error(
                "Internal node $(n.label) should not have an observed sequence: $(n.data)"
            )
        end
    end
    # nodes should not already have a reconstructed sequence
    for n in nodes(tree)
        if !isempty(n.data.state) || !isempty(n.data.reconstructed_sequence)
            error("Node $(n.label) has an already initialized fitch state: $(n.data)")
        end
    end
end

"""
    init_fitchstate!(dat::FitchData)
    init_fitchstate!(n::TreeNode)

Initialize the fitch state at `current_pos()` using the observed sequence.
If `dat.observed_sequence` is empty, the state remains empty.
If it is not empty but shorter than `current_pos()`, an error is thrown.
"""
function init_fitchstate!(dat::FitchData)
    if !isempty(dat.observed_sequence)
        if length(dat.observed_sequence) < current_pos()
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
    !isleaf(n) && pull_downstream_state!(n)
end
function pull_downstream_state!(n::TreeNode)
    # trying intersect
    # union with first child (to init.) then intersect with the rest
    c1, other_children = firstrest(children(n))
    union!(data(n).state, c1.data.state)
    for c in other_children
        intersect!(data(n).state, data(c).state)
    end
    # if intersect is empty, go for union!
    if isempty(data(n).state)
        for c in children(n)
            union!(data(n).state, data(c).state)
        end
    end

    return nothing
end

function fitch_down!(n::TreeNode)
    # get the state from ancestor
    if !isroot(n)
        pull_upstream_state!(n)
    end
    # sample state at current position
    sample_from_state!(n.data)
    # go to children
    foreach(fitch_down!, children(n))

    return nothing
end
function pull_upstream_state!(n::TreeNode)
    if length(ancestor(n).data.state) != 1
        throw(ErrorException("""
            Ancestor of $n does not have a fixed state on downward pass: $(a.data.state)
        """))
    end
    if isdisjoint(ancestor(n).data.state, n.data.state)
        # take the union
        union!(n.data.state, ancestor(n).data.state)
    else
        # intersection
        intersect!(n.data.state, ancestor(n).data.state)
    end
    return n.data.state
end

function sample_from_state!(dat::FitchData)
    if length(dat.reconstructed_sequence) != current_pos() - 1
        throw(ErrorException("""
            Trying to sample pos $(current_pos())\
            in sequence of length $(length(dat.reconstructed_sequence)).
            $dat
        """))
    end
    a = @chain length(dat.state) rand(1:_) nth(dat.state, _)
    push!(dat.reconstructed_sequence, a)
    dat.state = Set(a)
    return a
end



