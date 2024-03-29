"""
    fasta_to_tree(tree, fastafile)

Create a `Tree{FitchData}` object containing sequences of `fastafile` at leaves
"""
function fasta_to_tree(
    tree::Tree, fastafile::AbstractString;
)
    all_headers_in_tree = true
    all_leaves_in_fasta = true

    # get sequence type and copy tree with the appropriate data type
    tc = open(FASTA.Reader, fastafile) do reader
        data = first(reader) |> sequence |> FitchData
        tc = convert(Tree{typeof(data)}, tree)
    end

    #
    record = FASTA.Record()
    open(FASTA.Reader, fastafile) do reader
        while !eof(reader)
            read!(reader, record)
            if in(identifier(record), tc)
                tc[identifier(record)].data = FitchData(sequence(record))
            else
                all_headers_in_tree = false
            end
        end
    end

    for n in leaves(tc)
        if isempty(n.data.observed_sequence)
            all_leaves_in_fasta = false
            break
        end
    end
    !all_leaves_in_fasta && @warn "Not all leaves had a corresponding sequence \
    #     in the alignment (file: $fastafile)."
    !all_headers_in_tree && @warn "Some sequence headers in the alignment are \
    #     not found in the tree (file: $fastafile)."
    return tc
end

"""
    sequences_to_tree(tree::Tree, sequences::AbstractDict)

Create a `Tree{FitchData}` object containing sequences of `sequences` at leaves.
`sequences` should map a leaf label to a sequence.
"""
function sequences_to_tree(tree::Tree, sequences::AbstractDict)
    # creating tree with the right type
    tc = @chain sequences values first FitchData typeof convert(Tree{_}, tree)
    #
    for (label, seq) in Iterators.filter(x -> in(x[1], tc) && isleaf(tc[x[1]]), sequences)
        tc[label].data = FitchData(seq)
    end
    # check that all leaves have a sequence
    all_leaves_have_seq = all(n -> !isempty(n.data.observed_sequence), leaves(tc))
    !all_leaves_have_seq && @warn "Not all leaves had a corresponding sequence."
    return tc
end


