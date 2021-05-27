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

function evolve!(seq, w)
	for (i,nt) in enumerate(seq)
	 	newnt = sample(1:4, w)
	 	if newnt  != 1
	 		if nt == DNA_A
	 			if newnt == 2; seq[i] = DNA_C
	 			elseif newnt == 3; seq[i] = DNA_G
	 			elseif newnt == 4; seq[i] = DNA_T; end
	 		elseif nt == DNA_C
	 			if newnt == 2; seq[i] = DNA_G
	 			elseif newnt == 3; seq[i] = DNA_T
	 			elseif newnt == 4; seq[i] = DNA_A; end
	 		elseif nt == DNA_G
	 			if newnt == 2; seq[i] = DNA_T
	 			elseif newnt == 3; seq[i] = DNA_A
	 			elseif newnt == 4; seq[i] = DNA_C; end
	 		elseif nt == DNA_T
	 			if newnt == 2; seq[i] = DNA_A
	 			elseif newnt == 3; seq[i] = DNA_C
	 			elseif newnt == 4; seq[i] = DNA_G; end
	 		end
		end
	end
end
