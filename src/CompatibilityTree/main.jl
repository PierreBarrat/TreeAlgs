"""
	max_compatibility_splits(tab::CharTable)
	max_compatibility_splits(aln, labels)
	max_compatibility_splits(aln)

Splits of sequences that are compatible with the highest number of columns. A set of splits is compatible with a column of an alignment if the corresponding character can evolve without homoplasy (several mutations) on the tree defined by the splits.
"""
function max_compatibility_splits(tab::CharTable)
	if !isempty(tab.dat)
		cg = compatibility_graph(tab)
		clique = sort(LightGraphs.maximal_cliques(cg), by=x->length(x), rev=true)[1]
		return tree_popping(tab, clique)
	else
		return tree_popping(tab, 1:0)
	end
end
max_compatibility_splits(aln, labels) = max_compatibility_splits(CharTable(aln, labels))
max_compatibility_splits(aln) = max_compatibility_splits(CharTable(aln))

function tree_popping(tab::CharTable, sites)
	sp = sortperm(tab.labels)
	S = TreeTools.SplitList(tab.labels[sp])
	for i in sites
		s = get_split(tab.dat, i, sp)
		if !in(s, S.splits)
			push!(S.splits, s)
		end
	end
	sort!(S.splits, by=x->x.dat)
	return S
end

function get_split(X, i, perm)
	s = TreeTools.Split(size(X,1))
	@views for (m,x) in enumerate(X[perm,i])
		if x
			s.dat[m] = true
		end
	end
	return s
end
