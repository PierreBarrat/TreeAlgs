"""
	evolve!(t::Tree, L::Int, μ=1.; model = JC69(1.), seqkey=:seq)
"""
function evolve!(t::Tree, L::Int, μ=1.; model = JC69(1.), seqkey=:seq)
	t.root.data.dat[seqkey] = randseq(DNAAlphabet{4}(), L)
	for c in t.root.child
		evolve!(c, model, μ, seqkey)
	end
end

function evolve!(n::TreeNode, model, μ, seqkey=:seq)
	w = ProbabilityWeights(SubstitutionModels.P(model, μ * n.tau)[:,1], 1.)
	n.data.dat[seqkey] = deepcopy(n.anc.data.dat[seqkey])
	evolve!(n.data.dat[seqkey], w)
	for c in n.child
		evolve!(c, model, μ, seqkey)
	end
end

function evolve!(seq::NucSeq{4, DNAAlphabet{4}}, Q)
	Qt = Q'
	for (i, nt) in enumerate(seq)
		seq[i] = evolve(nt, Qt)
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
