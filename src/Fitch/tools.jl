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
    return nothing
end
