"""
	fasta_to_tree!(tree{MiscData}, fastafile::AbstractString, key=:seq; warn=true)

Add sequences of `fastafile` to nodes of `tree`.
For a leaf `n`, sequence is added to `n.data[key]`.
"""
function fasta_to_tree!(
	tree::Tree{MiscData}, fastafile::AbstractString, key = :seq;
	warn = true, default=missing,
)
	all_headers_in_tree = true
	all_leaves_in_fasta = true

	reader = open(FASTA.Reader, fastafile)
	record = FASTA.Record()
	while !eof(reader)
	    read!(reader, record)
	    if in(identifier(record), tree)
	    	tree.lnodes[identifier(record)].data[key] = sequence(record)
	    else
	    	all_headers_in_tree = false
	    end
	end
	close(reader)

	for n in leaves(tree)
		if !haskey(n.data, key)
			n.data[key] = default
			all_leaves_in_fasta = false
		end
	end
	warn && !all_leaves_in_fasta && @warn "Not all leaves had a corresponding sequence \
		in the alignment (file: $fastafile)."
	warn && !all_headers_in_tree && @warn "Some sequence headers in the alignment are \
		not found in the tree (file: $fastafile)."
	return nothing
end


"""
	write_fasta(file::AbstractString, tree, seqkey = :selfseq ; internal = false)
"""
function write_fasta(
	file::AbstractString, tree, seqkey = :seq;
	internal = false, root=false, sorted=true,
)
	open(FASTA.Writer, file) do f
		nodes_iter = sorted ? sort(collect(nodes(tree)), by = x->x.label) : nodes(tree)
		for n in nodes_iter
			if internal || (root && n.isroot) || n.isleaf
				rec = FASTA.Record(n.label, n.data[seqkey])
				write(f, rec)
			end
		end
	end
	return nothing
end
