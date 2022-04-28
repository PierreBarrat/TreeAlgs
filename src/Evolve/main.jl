"""
	evolve!(t::Tree, L::Int, μ=1.; model = JC69(1.), seqkey=:seq)

	initialize a random sequence of nucleotides of length `L` at the root, then simulate evolution
	by changing nucleotides at random according to a probability distribution specified 
	by the chosen `model`, which takes the mutation rate μ and branch lengths of nodes as parameters.

	The sequence at a node `n` can be accessed by `n.data.dat[seqkey]`, if a mutation has occured 
	on the branch between node `n` and its ancestor is given as a boolean in `n.data.dat["evolved"]`.
"""
function evolve!(t::Tree, L::Int, μ=1.; model = JC69(1.), seqkey=:seq)
	t.root.data.dat[seqkey] = randseq(DNAAlphabet{4}(), L)
	t.root.data.dat["evolved"]= false
	for c in t.root.child
		evolve!(c, model, μ, seqkey)
	end
end

function evolve!(n::TreeNode, model, μ, seqkey=:seq)
	n.data.dat["evolved"]= false
	w = ProbabilityWeights(SubstitutionModels.P(model, μ * n.tau)[:,1], 1.)
	n.data.dat[seqkey] = deepcopy(n.anc.data.dat[seqkey])
	evolve_seq!(n, w, seqkey)
	for c in n.child
		evolve!(c, model, μ, seqkey)
	end
end

function evolve_seq!(n::TreeNode, w, seqkey=:seq)
	seq = n.data.dat[seqkey]
	for (i,nt) in enumerate(seq)
	 	newnt = sample(1:4, w)
	 	if newnt  != 1
			n.data.dat["evolved"] = true
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
"""
	write_seq2fasta(t::Tree, fasta_name::String, output_dir::String, seqkey=:seq; only_terminals=false, remove_0_mutations=true)

	write evolved sequences to a fasta file with name `fasta_name`, `only_terminals` specifies if only terminal sequences should be 
	written to the file and `remove_0_mutations` if only sequences with mutations on their branches should be written to the file. 
"""
function write_seq2fasta(t::Tree, fasta_name::String, output_dir::String, seqkey=:seq; only_terminals=false, remove_0_mutations=true)
	mkpath(output_dir)
	open(FASTA.Writer, output_dir * "sequences" *fasta_name* ".fasta") do w
		if only_terminals
			iter = POTleaves(t)
		else
			iter = POT(t)
		end
		for node in iter
			if remove_0_mutations && !node.data.dat["evolved"]
				continue
			else
				x = node.data.dat[seqkey]
				rec = FASTA.Record(node.label, x)
				write(w, rec)
			end
		end
	end

end
