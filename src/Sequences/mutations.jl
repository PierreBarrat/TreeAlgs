"""
	struct Mutation{T}

```
	i::Int
	old::T
	new::T
```
"""
struct Mutation{T}
	i::Int
	old::T
	new::T
end
function Mutation(x::Tuple{Int,T,T} where T)
	return Mutation(x[1],x[2],x[3])
end

function Base.:(==)(x::Mutation, y::Mutation)
	mapreduce(f->getfield(x,f)==getfield(y,f), *, fieldnames(Mutation), init=true)
end
Base.hash(x::Mutation, h::UInt) = hash((x.i, x.old, x.new), h)

Base.reverse(x::Mutation) = Mutation(x.i, x.new. x.old)
function isreverse(x::Mutation, y::Mutation)
	(x.i == y.i) && (x.old == y.new) && (x.new == y.old)
end


function parse_mutation(mut::AbstractString)
	oldstate = mut[1]
	newstate = mut[end]
	pos = parse(Int, mut[2:end-1])
	return Mutation(pos, oldstate, newstate)
end
function parse_mutation(mut::AbstractString, ::Val{T}) where T
	oldstate = mut[1]
	newstate = mut[end]
	pos = parse(Int, mut[2:end-1])
	return Mutation{T}(pos, oldstate, newstate)
end
parse_mutation(mut::AbstractString, T) = parse_mutation(mut, Val(T))
parse_mutation(mut::Mutation) = mut


"""
	mutations_from_sequences!(t::Tree{MiscData}; seqkey=:seq, mutkey=:muts)
"""
function mutations_from_sequences!(t::Tree{MiscData}; seqkey=:seq, mutkey=:muts)
	t.root.data[mutkey] = []
	for n in Iterators.filter(!isroot, nodes(t))
		n.data[mutkey] = []
		for (i, a) in enumerate(n.data[seqkey])
			if a != n.anc.data[seqkey][i]
				push!(n.data[mutkey], Mutation(i, n.anc.data[seqkey][i], a))
			end
		end
	end

	return nothing
end
