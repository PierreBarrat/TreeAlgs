using FASTX
using TreeTools

"""
	fasta2tree!(tree, fastafile::AbstractString, key=:seq; warn=true)

Add sequences of `fastafile` to nodes of `tree`.
For a leaf `n`, sequence is added to `n.data.dat[key]`.
If `key` is a `Tuple`, store sequence in `n.data.dat[key[1]][key[2]]...
"""
function fasta2tree!(tree, fastafile::AbstractString, key::Union{Symbol, AbstractString}=:seq; warn=true)
	reader = open(FASTA.Reader, fastafile)
	record = FASTA.Record()
	while !eof(reader)
	    read!(reader, record)
	    if haskey(tree.lnodes, identifier(record))
	    	tree.lnodes[identifier(record)].data.dat[key] = sequence(record)
	    end
	end
	#
	flag = true
	for (name, n) in tree.lleaves
		if !haskey(n.data.dat, key)
			n.data.dat[key] = missing
			flag = false
		end
	end
	warn && !flag && @warn "Not all leaves had a corresponding sequence in the alignment (file: $fastafile)."
	return flag
end


function fasta2tree!(tree, fastafile::AbstractString, ks::Tuple; warn=true)
	# Setting dicts
	for n in values(tree.lleaves)
		recursive_key_init!(n.data.dat, ks[1:end-1]...)
	end
	key = ks[end]
	#
	reader = open(FASTA.Reader, fastafile)
	record = FASTA.Record()
	while !eof(reader)
	    read!(reader, record)
	    if haskey(tree.lnodes, identifier(record))
	    	recursive_set!(tree.lnodes[identifier(record)].data.dat, sequence(record), ks...)
	    end
	end
	#
	flag = true
	for (name, n) in tree.lleaves
		if !recursive_haskey(n.data.dat, ks...)
			recursive_set!(n.data.dat, missing, ks...)
			flag = false
		end
	end
	warn && !flag && @warn "Not all leaves had a corresponding sequence in the alignment (file: $fastafile)."
	return flag
end




"""
	write_fasta(file::AbstractString, tree, seqkey = :selfseq ; internal = false)
"""
function write_fasta(
	file::AbstractString, tree, seqkey = :selfseq;
	internal = false, root=false, sorted=true,
)
	open(FASTA.Writer, file) do f
		nodes_iter = sorted ? sort(collect(nodes(tree)), by = x->x.label) : nodes(tree)
		for n in nodes_iter
			if internal || (root && n.isroot) || n.isleaf
				rec = FASTA.Record(n.label, recursive_get(n.data.dat, seqkey))
				write(f, rec)
			end
		end
	end
	return nothing
end

##################################################
##### Utilities for dicts with nested keys #######
##################################################
function recursive_key_init!(dat, key, ks...)
    if !haskey(dat, key)
        dat[key] = Dict()
    end
    recursive_key_init!(dat[key], ks...)
end
recursive_key_init!(dat) = nothing
function recursive_get(dat, key, ks...)
    if isempty(ks)
        return dat[key]
    else
        return recursive_get(dat[key], ks...)
    end
end
recursive_get(dat, key::Tuple) = recursive_get(dat, key...)
function recursive_set!(dat, value, key, ks...)
    if isempty(ks)
        dat[key] = value
    else
        if !haskey(dat, key)
            dat[key] = Dict()
        end
        recursive_set!(dat[key], value, ks...)
    end
    dat
end
recursive_set!(dat, value, key::Tuple) = recursive_set!(dat, value, key...)
function recursive_push!(dat, value, key, ks...)
    if isempty(ks)
        push!(dat[key], value)
    else
        recursive_push!(dat[key], value, ks...)
    end
    dat
end
recursive_push!(dat, value, key::Tuple) = recursive_push!(dat, value, key...)
function recursive_haskey(dat, key, ks...)
    if !haskey(dat, key)
        return false
    elseif isempty(ks)
        return true
    else
        return recursive_haskey(dat[key], ks...)
    end
end

