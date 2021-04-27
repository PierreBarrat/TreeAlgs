struct FitchState{T}
	state::Set{T}
end
FitchState(::Val{T}) where T  = FitchState{T}(Set{T}())
FitchState(a::T) where T = FitchState{T}(Set(a))
isempty(fs::FitchState) = isempty(fs.state)
length(s::FitchState) = length(s.state)


"""
	fitch!(t::Tree, outkey=:ancestral_seq, seqkey=:seq; clear_fitch_states=true, variable_positions=missing)
"""
function fitch!(t::Tree, outkey=:seq, seqkey=:seq; clear_fitch_states=true, variable_positions=Int64[])
	fitchkey = :fitchstate
	# Initializing ancestral sequences
	seq = TreeTools.recursive_get(first(t.lleaves)[2].data.dat, seqkey)
	init_ancestral_sequences!(t, outkey, seq)

	#Getting variable positions
	if ismissing(variable_positions)
		variable_positions = 1:length(seq)
	elseif isempty(variable_positions)
		variable_positions = get_variable_positions(t, seqkey)
	end

	# Algorithm
	for i in 1:length(seq)
		if in(i, variable_positions)
			init_fitchstates!(t, i, seqkey, fitchkey)
			fitch_up!(t, fitchkey)
			fitch_remove_gaps!(t, fitchkey)
			fitch_root_state!(t, fitchkey)
			fitch_down!(t, fitchkey)
			fitch_remove_gaps!(t, fitchkey)
			fitch_sample!(t, outkey, fitchkey)
		else
			for n in values(t.lnodes)
				if !n.isleaf
					TreeTools.recursive_push!(n.data.dat, seq[i], outkey)
				end
			end
		end
	end
	# Clearing fitch states
	if clear_fitch_states
		for n in values(t.lnodes)
			delete!(n.data.dat, :fitchstate)
		end
	end
end

function init_ancestral_sequences!(t, outkey, seq)
	for n in values(t.lnodes)
		if !n.isleaf
			n.data.dat[outkey] = similar(seq, 0)
		end
	end
end
function init_ancestral_sequences!(t, outkey::Tuple, seq)
	for n in values(t.lnodes)
		if !n.isleaf
			TreeTools.recursive_key_init!(n.data.dat, outkey[1:end-1]...)
			TreeTools.recursive_set!(n.data.dat, similar(seq, 0), outkey...)
		end
	end
end

"""
	init_fitchstates!(t::Tree, i::Int64, seqkey = :seq, fitchkey=:fitchstate)
"""
function init_fitchstates!(t::Tree, i::Int64, seqkey::Union{Symbol, AbstractString}=:seq, fitchkey=:fitchstate)
	for n in values(t.lleaves)
		n.data.dat[fitchkey] = FitchState(n.data.dat[seqkey][i])
	end
end
function init_fitchstates!(t::Tree, i::Int64, seqkey::Tuple, fitchkey=:fitchstate)
	for n in values(t.lleaves)
		n.data.dat[fitchkey] = FitchState(TreeTools.recursive_get(n.data.dat, seqkey...)[i])
	end
end

"""
	ancestral_state(fstates::Vararg{FitchState{T}}) where T
"""
function ancestral_state(fs::FitchState{T}, fstates::Vararg{FitchState{T}}) where T
	aFs = FitchState(Val(T))
	union!(aFs, fs)
	for s in fstates
		intersect!(aFs, s)
	end
	if isempty(aFs) || (length(aFs.state) == 1 && isgap(first(aFs.state)))
		union!(aFs, fs)
		for s in fstates
			union!(aFs, s)
		end
	end
	return aFs
end

"""
	get_downstream_state!(an::TreeNode, fitchkey=:fitchstate)
"""
function get_downstream_state!(an::TreeNode, fitchkey=:fitchstate)
	an.data.dat[fitchkey] = ancestral_state((n.data.dat[fitchkey] for n in an.child)...)
	nothing
end

"""
"""
function fitch_up!(r::TreeNode, fitchkey)
	for c in r.child
		!c.isleaf && fitch_up!(c, fitchkey)
	end
	get_downstream_state!(r, fitchkey)
end
fitch_up!(t::Tree, fitchkey=:fitchstate) = fitch_up!(t.root, fitchkey)

"""
"""
function fitch_root_state!(t::Tree, fitchkey=:fitchstate)
	fs = t.root.data.dat[fitchkey]
	L = Dict{Any,Float64}()
	# Compute likelihood of each possible state 
	for (k,a) in enumerate(fs.state) 
		for c in t.root.child
			if !ismissing(c.data.tau)
				!haskey(L,a) && (L[a] = 0.)
				if isempty(intersect([a], c.data.dat[fitchkey].state))
					# Mutation needed
					L[a] += 1 - exp(-c.data.tau)
				else
					L[a] += exp(-c.data.tau)
				end
			end
		end
	end
	amax = isempty(L) ? rand(fs.state) : findmax(L)[2]
	filter!(==(amax), fs.state)
end

"""
	set_child_state!(fs_child::FitchState{T}, fs_anc::FitchState{T}) where T
"""
function set_child_state!(fs_child::FitchState{T}, fs_anc::FitchState{T}) where T
	if length(fs_anc) > 1
		error("State $(fs_anc.state) on downward pass")
	end
	if in(first(fs_anc.state), fs_child.state)
		filter!(==(first(fs_anc.state)), fs_child.state)
	else
		a = rand(fs_child.state)
		filter!(==(a), fs_child.state)
	end
end
"""
"""
function get_upstream_state!(n::TreeNode, fitchkey=:fitchstate)
	set_child_state!(n.data.dat[fitchkey], n.anc.data.dat[fitchkey])
end
function fitch_down!(r::TreeNode, fitchkey)
	!r.isroot && get_upstream_state!(r, fitchkey)
	for c in r.child
		!c.isleaf && fitch_down!(c, fitchkey)
	end
end
fitch_down!(t::Tree, fitchkey=:fitchstate) = fitch_down!(t.root, fitchkey)

"""
"""
function fitch_remove_gaps!(t, fitchkey=:fitchstate)
	for n in values(t.lnodes)
		if !n.isleaf && length(n.data.dat[fitchkey]) > 1
			filter!(!isgap, n.data.dat[fitchkey].state)
		end
	end
end

"""
"""
function fitch_sample!(t::Tree, outkey::Tuple, fitchkey=:fitchstate)
	for n in values(t.lnodes)
		if !n.isleaf
			TreeTools.recursive_push!(n.data.dat, rand(n.data.dat[fitchkey].state), outkey...)
		end
	end
end
function fitch_sample!(t::Tree, outkey::Union{Symbol, AbstractString}, fitchkey=:fitchstate)
	for n in values(t.lnodes)
		if !n.isleaf
			push!(n.data.dat[outkey], rand(n.data.dat[fitchkey].state))
		end
	end
end
fitch_sample(fs::FitchState{DNA}) = LongDNASeq([rand(s) for s in fs.state])
fitch_sample(fs::FitchState{RNA}) = LongRNASeq([rand(s) for s in fs.state])
fitch_sample(fs::FitchState{AminoAcid}) = LongAminoAcidSeq([rand(s) for s in fs.state])
fitch_sample(fs::FitchState{Char}) = LongCharSeq([rand(s) for s in fs.state])