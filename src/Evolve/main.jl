"""
	evolve!(t::Tree{TreeTools.MiscData}, L::Int, μ=1.; model = JC69(1.), seqkey=:seq)

Evolve DNA sequences of length `L` along the branches of `t` using `model`.
The sequence at the root is random: `randseq(DNAAlphabet{4}(), L)`.
For each node `n` of `t`, the sequence is stored in `n.data[seqkey]`.
"""
function evolve!(t::Tree{TreeTools.MiscData}, L::Int, μ=1.; model = JC69(1.), seqkey=:seq)
	t.root.data[seqkey] = randseq(DNAAlphabet{4}(), L)
	t.root.data.dat["evolved"] = false
	for c in t.root.child
		evolve!(c, model, μ, seqkey)
	end
end

function evolve!(n::TreeNode{TreeTools.MiscData}, model, μ::Real, seqkey=:seq)
	w = ProbabilityWeights(SubstitutionModels.P(model, μ * branch_length(n))[:,1], 1.)
	n.data[seqkey] = deepcopy(n.anc.data.dat[seqkey])
	evolve!(n.data.dat[seqkey], w)
	for c in n.child
		evolve!(c, model, μ, seqkey)
	end
end

function evolve!(n::TreeNode{TreeTools.MiscData}, Q)
	Qt = Q'
	seq = n.data.dat[seqkey]
	n.data.dat["evolved"]= false
	for (i, nt) in enumerate(seq)
		seq[i] = evolve(nt, Qt)
		seq[i]=!nt && (n.data.dat["evolved"]=true) ##check if sequence has changed
	end
	return nothing
end

function evolve(nt::DNA, Qt)
	nts = [DNA_A, DNA_C, DNA_G, DNA_T]
	if !in(nt, nts)
		return nt
	else
		s = onehot(nt)
		w = Qt * s
		return sample(nts, ProbabilityWeights(w))
	end
end

function onehot(c::DNA)
	s = Int8[0, 0, 0, 0]
	if c == DNA_A
		s[1] = 1
	elseif c == DNA_C
		s[2] = 1
	elseif c == DNA_G
		s[3] = 1
	elseif c == DNA_T
		s[4] = 1
	end
	return s
end
